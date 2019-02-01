.libPaths(c(.libPaths(), "/mnt/work1/users/pughlab/bin/r_lib/library/3.1"))
###############################################################
#
#                 CNV Fingerprinter: All-vs-All
#
# Author: Rene Quevedo
# Last Modified: July-2-2015
# Date Created: March-05-2015
###############################################################
# Function: Creates and plots an all-against-all matrix of 
#   cell lines from all datasets against all cell lines of the
#   same dataset and all others as well.  
#
#   Creates a heatmap of percent identity based on SSE by taking
#   the ratio of binned segments between two files and seeing how
#   much that deviates from a perfect 1 ratio.
###############################################################
# Arguments:
#   [1]  Output Directory
#   [2]  Location of all the input ABSOLUTE SEG_MAF files
#   [3]  Reading the segtab files in for the first time (Default = 1)
#   [4]  Generate all individual comparison plots (Default = 1)
#   [5]  Location of the ucsc.chromInfo file
#   [6]  Location containing the cell line annotation file   
#   [7]  Chr Segment Size  (Default = 20000)
#   [8]  Plot List - Give the segtab files to plot specifically if not all output is wanted
###############################################################
# Example Execution:
#   Rscript cnvAvaFingerprinter.R \
#[1]     /mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/hapmap_norm/cnv_fingerprint/fingerprint/output \
#[2]     /mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/hapmap_norm/cnv_fingerprint/absolute_output/total/SEG_MAF \
#[3]     1 \
#[4]     1 \
#[5]     /mnt/work1/users/pughlab/references/ucsc.hg19.chromInfo.txt \
#[6]     /mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/hapmap_norm/annotations/merged_annotations.Rdata \  
#[7]     5000
###############################################################
library(gplots)     # Creating ggplots2 plots
library(gtools)     # Used to "smartbind" the dataframes together
library(plyr)       # For counting high SD bins
library(doParallel) # For Parallelization
library(foreach)    # Parallelization as well
library(intervals)  # Dealing with Intervals and intersecting intervals (CN segments)
library(data.table) # For fast-order joining of many matrices and lists
library(pastecs)    # Used to convert the density plots into time-series
library(corrgram)   # For creating Correlograms
library(scales)     # For identifying the turnpoints in a density plot (local minimas, maximas) and alpha colouring
library(doSNOW)


#########################################################################################################
#         Variables
#########################################################################################################
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)


cn.dir <- '/Users/rquevedo/Desktop/bhk_lab/results/cnv_fingerprinting/test.subset/hscr.a2/segtab'
#ucsc.chromInfo <- '/Users/rquevedo/Desktop/bhk_lab/reference/ucsc.hg19.chromInfo.txt'
ucsc.chromInfo <- '/mnt/work1/users/pughlab/references/ucsc.hg19.chromInfo.txt'
cl.anno.file <- "/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/hapmap_norm/annotations/merged_annotations.Rdata"
#cl.anno.file <- "/Users/rquevedo/Desktop/bhk_lab/annotations/merged_annotations.Rdata"
global.chrom <- c(1:22)
global.bin.size <- 50000 #http://dgv.tcag.ca/dgv/app/statistics?ref=
max.diff <- 10   #Used to calculate the normalization coefficient - average max difference between two profiles
mismatch.val <- 0

# Shortcut the Data Flags
generate.output <- 1  # Flag to generate the output files
process.data <- 1 # Flag to process the data
ignore.process <- 0  # Flag to ignore processing and just print out individual plots

#Reading in the files and pre-processing them:
output.dir <- args[1]
cn.dir <- args[2]
rdata.dir <- "rdata"
create.df <- 1  # 1 = Read each cn.file for the first time
                # 0 = Use pre-existing cn.file df Rdata files
existing.rdata.dir <- '';  # Set the directory if create.df is set to 0, else default "./rdata" will be used
if(length(args) >= 3){
  create.df <- as.numeric(args[3])
}

# Generating the CN comparison plots:
plots.dir <- "plots"
gen.plots <- 1  # 1 = Generate plots
                # 0 = Don't generate plots
plot.scatterplot <- 0  # 1 = Generates the scatterplot and ellipse correlogram  (Warning: VERY lengthy procedure)
                       # 0 = Don't generate the correlogram plot 
if(length(args) >= 4){
  gen.plots <- as.numeric(args[4])
}
min.thresh <- 1.5
max.thresh <- 3

# Extra optional arguments:
if(length(args == 7)){
  ucsc.chromInfo <- args[5]
  cl.anno.file <- args[6]
  global.bin.size <- as.numeric(args[7])
}
plot.list <- c()
if(length(args) == 8){
  plot.list.file <- args[8]
  plot.list.df <- read.csv(plot.list.file, header = F)
  plot.list <- as.vector(plot.list.df[,1])
}

#Parallel Processing
num.of.cores <- 4 #Samwise has 40
dec.cores <- 2    # if limit is hit, decrease number of required cores by this amount

# Headers to scan for:
# Copy ratio values used to compare
cr.header1 <- 'hscr.a1'
cr.header2 <- 'hscr.a2'
header.list <- list(cr.header1, cr.header2)
# Chrom location of segment
chrom.header <- 'Chromosome'
start.header <- 'Start.bp'
end.header <- 'End.bp'
# Annotation File header
fn.headers <- c('su2c.filename', 'cgp.filename', 'ccle.filename', "pfizer.filename.x", "pfizer.filename.y", 'pfizer.filename', 'gdsc.filename.x', 'gdsc.filename.y')

# Output Files:
out.cn.sse <- "cn_sse.Rdata"
out.cn.files <- "row_col_filenames.Rdata"
out.total.rdata <- "total_rdata_list.Rdata"
out.spearman.mat <- "spearman_matrix.Rdata"
out.eucl.mat <- "euclidean_matrix.Rdata"


print(paste(Sys.time(), ": Reading in UCSC Chrom Info...", sep=""))
ucsc.chrom.df <- read.csv(ucsc.chromInfo, header=T, sep="\t", check.names=FALSE)
ucsc.chrom.df$chrom <- gsub("^chr", "", ucsc.chrom.df$chrom)



#########################################################################################################
#         Functions
#########################################################################################################
# Function: compareCnProfiles
# Purpose:  Uses the merged annotation file to retrieve all matching names for every dataset
# Input:  The name of the input cel file  
# Returns:  List: Contains all names of matching cel files for all datasets.
compareCnProfiles <- function(cn.binned.df1, cn.binned.df2, seg.col.name, num.bins, mdiff, 
                              min.thresh=1.5, max.thresh=3){
  cn1 <- cn.binned.df1[,seg.col.name]
  cn2 <- cn.binned.df2[,seg.col.name]
  
  #SSD Value (in -log space, normalized to worst case scenario)
  cn.squared.diff <- (abs(cn1 - cn2)+1)^2
  normalization.coefficient <- num.bins * mdiff
  cn.ssd <- -log(sum(cn.squared.diff) / normalization.coefficient)
  #cn.ssd <- sum(cn.squared.diff)
  
  #Density function to find first local minima
  dens.cn.sd <- density(cn.squared.diff)
  ts_y<-ts(dens.cn.sd$y)
  #Identify turnpoints using the pastecs package
  tp<-turnpoints(ts_y)
  loc.min <- dens.cn.sd$x[tp$tppos[2]]
  if(is.na(loc.min)){
    loc.min <- 1
  } else if(loc.min < min.thresh | loc.min > max.thresh){
    loc.min <- min.thresh
  } 
  
  
  # index which bars are coloured for loss, gain, or non-significant
  thresh.cn.sd <- cn.squared.diff[which(cn.squared.diff > loc.min)]
  thresh.cn.ssd <- -log(sum(thresh.cn.sd) / normalization.coefficient)
  
  cn.ssd.list <- list(cn.sd = cn.squared.diff, 
                      ssd = cn.ssd, 
                      thresh.cn.sd = thresh.cn.sd, 
                      thresh.ssd = thresh.cn.ssd, 
                      loc.min = dens.cn.sd$x[tp$tppos[2]])  
  return(cn.ssd.list)
}


getSsdScoreTally <- function(cn.binned.df1, cn.binned.df2, seg.col.name, num.bins, mdiff){
  cn1 <- cn.binned.df1[,seg.col.name]
  cn2 <- cn.binned.df2[,seg.col.name]
  
  
  #SSD Value (in -log space, normalized to worst case scenario)
  cn.squared.diff <- (abs(cn1 - cn2)+1)^2
  normalization.coefficient <- num.bins * mdiff
  cn.ssd <- -log(sum(cn.squared.diff) / normalization.coefficient)
  #cn.ssd <- sum(cn.squared.diff)
  
  return(list(ssd = cn.ssd))
}

#Selectively replaces "CEL$" (case insensitive) with "[seg.header].Rdata"
formatFileList <- function(file.l, seg.header){
  replace.vector <- grep("cel$", file.l, ignore.case=TRUE)
  file.l[replace.vector] <- gsub("cel$", paste(seg.header, ".Rdata", sep=""), ignore.case = TRUE, file.l[replace.vector])
  return(file.l)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
  
fillSimilarityMatrix <- function(sim.matrix){
  for(each.col in c(1:dim(sim.matrix)[1])){
    for(each.row in c(1:dim(sim.matrix)[2])){
      if(is.na(sim.matrix[each.row, each.col])){
        val <- sim.matrix[each.col, each.row] 
        sim.matrix[each.row, each.col] <- val
      }
    }
  }
  return(sim.matrix)
}

plotSd <- function(file.l=NULL, comp.all=FALSE, total.rdata.list=total.rdata.list, seg.header=NULL, plots.path=NULL){
  require(pastecs)    
  require(corrgram)
  require(scales)
  
  comp.against <- file.l
  
  #Ensure proper format of arguments to the subfunction
  if(is.null(file.l)){
    stop("Error: A file list was not provided for SD-plotting")
  } else if (is.null(seg.header)){
    stop("Error: No seg header was provided")
  } else if (is.null(plots.path)){
    stop("Error: No output directory for plots was given")
  }
  if(length(file.l[which(!file.l %in% names(total.rdata.list))]) > 0){
    warning(paste("File(s): ", file.l[which(!file.l %in% names(total.rdata.list))], 
                  "were not found in the total.rdata.list", sep=""))
  }
  if(comp.all==TRUE){
    print("Comparing given files to all files in total.rdata.list")
    comp.against <- names(total.rdata.list)
  }
  
  # Removes any possible .CEL(case-insensitive) heading and replaces it with seg.header
  file.l <- formatFileList(file.l, seg.header)
  
  # Cycles through all comparing files to generate plots
  for(file.a in file.l){
    file.a.index <- which(names(total.rdata.list) %in% file.a)
    for(file.b in comp.against){
      file.b.index <- which(names(total.rdata.list) %in% file.b)
      
      if(!file.a == file.b){
        print("Comparing CN profiles to generate SD plots")
        ab.sd.list <- compareCnProfiles(total.rdata.list[[file.a.index]], total.rdata.list[[file.b.index]], 
                                        seg.header, total.bins, max.diff, min.thresh=2)
        
        #setting up some of the plotting parameters
        cn1 <- total.rdata.list[[file.a.index]][,seg.header]
        cn2 <- total.rdata.list[[file.b.index]][,seg.header]
        max.ylim.sd <- 15
        max.ylim <- 3
        
        dir.create(file.path(plots.path, file.a), showWarnings=F, recursive=T) 
        pdf(file.path(plots.path, file.a, gsub(".Rdata", ".pdf", file.b)))
        print("Generating 'Plot 1': Comparative HSCR segmented chromosomal plot...")
        plotSegmentedHscr(cn1=cn1, cn2=cn2, max.ylim=max.ylim, max.ylim.sd=max.ylim.sd, fname1=file.a, fname2=file.b)
        
        diff.cn <- cn2-cn1
        print("Generating 'Plot 2': Differences and thresholds between CN profile 1 and 2...")
        plotSegmentedSd(comp.list=ab.sd.list, diff.cn=diff.cn, 
                        max.ylim=max.ylim, max.ylim.sd=max.ylim.sd)
        dev.off()
      }
      
    }
  }
}


plotSegmentedHscr <- function(cn1=NULL, cn2=NULL, max.ylim=NULL, max.ylim.sd=NULL, fname1=NULL, fname2=NULL){
  par(mfrow=c(1,1))
  
  #Creates a blank plot
  plot(c(1:length(cn1)), 
       type="n", axes=T, 
       ylim=c(0,max.ylim), 
       main=fname1, cex.main=0.7,
       xlab="Genomic Bins", cex.lab=0.7,
       ylab="HSCR",
       cex.axis=0.7)
  gen.pos.match <- which(cn1==cn2)
  hscr.match <- cn1[cn1==cn2]
  # Adds individual and matching dots
  points(c(1:length(cn1)), cn1, pch=20, col=alpha("black", 0.5))
  points(c(1:length(cn2)), cn2, pch=20, col=alpha("darkorange", 0.5))
  points(gen.pos.match, hscr.match, pch=20, col=alpha("gold2", 0.5))
  
  addChrRect(max.ylim.sd=max.ylim.sd)
  
  # Creates a legend
  legend(30000,3,
         c(fname1,fname2),
         lty=c(1,1), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5),col=c("black", "darkorange"), # gives the legend lines the correct color and width
         cex=0.40)
}


plotSegmentedSd <- function(comp.list=NULL, diff.cn=NULL, max.ylim=NULL, max.ylim.sd=NULL){
  #### PLOT #3: #Squared-Diff plot between CN profile 1 and 2
  bar.cols <- c("blue", "red", "black")
  layout(matrix(c(1,1,1,1,0,2,2,2,2,3), 2, 5, byrow = TRUE))
  
  # index which bars are coloured for loss, gain, or non-significant
  pos <- diff.cn
  pos[which(comp.list$cn.sd <= comp.list$loc.min)] <- 0
  pos[which(pos > 0)] <- 2
  pos[which(pos < 0)] <- 1
  pos[which(pos == 0)] <- 3
  
  # Plot 3a) Difference from 2nd CN bins to the reference
  par(mar=c(0,4,2,0))
  barplot(diff.cn,
          col = bar.cols[pos], border = bar.cols[pos],
          ylab="Difference in CN-segments",
          xlab="Genomic bins")
  addChrRect('barplot', max.ylim.sd=max.ylim.sd)
  
  
  # Plot 3b) Cubed difference barplots
  par(mar=c(4,4,1,0))
  barplot(comp.list$cn.sd,
          col = bar.cols[pos], border = bar.cols[pos],
          ylim=c(0,max.ylim.sd),
          xlab="Genomic Bins",
          ylab="Squared-Difference in CN-segments")
  abline(comp.list$loc.min, b=0, col="red")
  addChrRect('barplot', max.ylim.sd=max.ylim.sd)
  
  
  dens.cn.sd <- density(comp.list$cn.sd)
  # Plot 3c) Density plot of cubed CN segments to identify noise
  par(mar=c(3.4,0,1,2))
  plot(dens.cn.sd$y, dens.cn.sd$x, 
       type="n", axes=T,
       xlab='', ylab='',
       ylim=c(0,max.ylim.sd),
       yaxt='n', xaxt='n')
  lines(dens.cn.sd$y, dens.cn.sd$x)
  abline(comp.list$loc.min, b=0, col="red")
  layout(c(1))
}


# Function: populateGenomeBins
# Purpose:  Takes in a pre-determined blank list of all bin sizes for each chromosome in the human genome
#   and maps out the given segment length and CN-segments to each of the bins.  This will allow for an easier
#   bin-vs-bin comparison in future steps.
# Input:  cn.val.f1: dataframe containing the CN segments and their seg-values (hscr.a1, seg_mean, etc..)
#         chr.bins: The blank list containing all preformatted empty bins
#         seg.col.name: column name for the segment header (i.e. seg_mean, hscr.a1, hscr.a2, etc...)
# Returns:  genome.tbins.df: Dataframe containing all the preformatted bins with their associated seg-values 
#   according to cn.val.f1.
populateGenomeBins <- function(cn.val.f1, chr.bins, seg.col.name){
  genome.tbins.df <- data.frame(matrix(ncol=5, nrow=0))
  for(each.chr in global.chrom){
    composite.bin.df <- as.data.frame(matrix(ncol=6, nrow=0)) #Keeps track of segments that do not fit whole bins
    matching.intervals.df <- as.data.frame(matrix(ncol=4, nrow=0)) # Keeps track of whole bin segments.
    
    colnames(matching.intervals.df) <- c("Start.bin", "End.bin", seg.col.name, 'copy.ratio')
    blank.chr.int <- Intervals(chr.bins[[each.chr]], closed=TRUE)
    cn.chr.f1 <- cn.val.f1[cn.val.f1$Chromosome %in% each.chr,]
    for(each.chr.row in 1:dim(cn.chr.f1)[1]){
      cn.int.f1 <- Intervals(matrix(c(cn.chr.f1[each.chr.row, 'Start.bp'], cn.chr.f1[each.chr.row, 'End.bp']), ncol=2))
      match.matrix <- as.matrix(interval_intersection(cn.int.f1, blank.chr.int))
      # Matches the blank matrix to the "match matrix" that has removed first and last row (to be calculated later)
      temp.matrix <- matrix(sort(intersect(chr.bins[[each.chr]], match.matrix[-c(1,nrow(match.matrix)),])), ncol=2, byrow=T)
      temp.matrix <- cbind(temp.matrix, rep(cn.chr.f1[each.chr.row, seg.col.name], nrow(temp.matrix)))
      temp.matrix <- cbind(temp.matrix, rep(cn.chr.f1[each.chr.row, 'copy.ratio'], nrow(temp.matrix)))
      colnames(temp.matrix) <- c("Start.bin", "End.bin", seg.col.name, 'copy.ratio')
      matching.intervals.df <- rbind(matching.intervals.df, temp.matrix)
      colnames(matching.intervals.df) <- c("Start.bin", "End.bin", seg.col.name, 'copy.ratio')
      
      # Stores the intervals that do not fill a complete bin in the composite dataframe to handle after.
      # bp-level intervals match in the bin
      first.match.row <- match.matrix[1,]
      last.match.row <- match.matrix[nrow(match.matrix),]
      # the bin being matched to
      first.match.blank.int <- Intervals(first.match.row)
      last.match.blank.int <- Intervals(last.match.row)
      composite.matrix <- matrix(c(first.match.row, blank.chr.int[!is.na(interval_overlap(blank.chr.int, first.match.blank.int) > 0)],
                                   last.match.row, blank.chr.int[!is.na(interval_overlap(blank.chr.int, last.match.blank.int) > 0)]), 
                                 nrow=2, byrow=T)
      composite.matrix <- cbind(composite.matrix, 
                                rep(cn.chr.f1[each.chr.row, seg.col.name], nrow(composite.matrix)),
                                rep(cn.chr.f1[each.chr.row, 'copy.ratio'], nrow(composite.matrix)))
      if(length(which(first.match.row == last.match.row)) ==  2){   # If the first.match.row vector matches last.match.row vector
        composite.bin.df <- data.frame(rbind(as.matrix(composite.bin.df), composite.matrix[1,]))
      } else {
        composite.bin.df <- data.frame(rbind(as.matrix(composite.bin.df), composite.matrix))
      }
      # Composite matrix result: start.bp, end.bp, start.bin, end.bin, segment.col.name
      
    }
    colnames(composite.bin.df) <- c("Start.bp", "End.bp", "Start.bin", "End.bin", seg.col.name, 'copy.ratio')
    matching.intervals.df <- rbind(matching.intervals.df, generateWeightedCompositeBin(composite.bin.df, seg.col.name))
    matching.intervals.df <- matching.intervals.df[order(matching.intervals.df$Start.bin),]
    #Finds the intervals not included in the analysis and fills them in
    match.intervals <- Intervals(matrix(c(matching.intervals.df[, 'Start.bin'], 
                                          matching.intervals.df[, 'End.bin']), ncol=2, byrow=F))
    mismatched.int.matrix <- as.matrix(interval_difference(blank.chr.int, match.intervals))
    mismatched.int.matrix<- cbind(mismatched.int.matrix, 
                                  rep(mismatch.val, dim(mismatched.int.matrix)[1]), 
                                  rep(mismatch.val, dim(mismatched.int.matrix)[1]))
    matching.intervals.df <- data.frame(rbind(as.matrix(matching.intervals.df), mismatched.int.matrix))
    
    matching.intervals.df <- matching.intervals.df[order(matching.intervals.df$Start.bin),]
    matching.intervals.df$chr <- each.chr
    
    
    
    genome.tbins.df <- rbind (genome.tbins.df, matching.intervals.df)
    # Need to find the complement of the intervals using interval_complement
  }
  colnames(genome.tbins.df) <- c("Start.bin", "End.bin", seg.col.name, 'copy.ratio', "Chr")
  return(genome.tbins.df)
}

# Function: generateWeightedCompositeBin
# Purpose:  Looks at predetermined bins and proportionally weighs each CN-segment that makes up that bin
#   to return a weighted-CN-segment for that bin.
# Input:  composite.bin.df <- dataframe containing the actual intervals that make up each predetermined bins
#         seg.col.name <- column name for the segment header (i.e. seg_mean, hscr.a1, hscr.a2, etc...)
# Returns:  composite.bin.final.df: data-frame containing the genome binned into the set sizes and having a 
#   proportionally weighted CN segment according to to the segment length in that bin.
generateWeightedCompositeBin <- function(composite.bin.df, seg.col.name){
  composite.bin.final.df <- as.data.frame(matrix(ncol=3, nrow=0))
  for(each.uniq.bin in unique(composite.bin.df$Start.bin)){
    unique.bin.df <- composite.bin.df[composite.bin.df$Start.bin == each.uniq.bin, ]
    
    start.bin <- each.uniq.bin
    end.bin <- unique.bin.df[1,'End.bin']
    
    # Calculates the proportionally-weighted CN segment for the given bins    
    unique.bin.df$length <- unique.bin.df$End.bp - unique.bin.df$Start.bp
    total.cov <- sum(unique.bin.df$length)
    unique.bin.df$weighted.seg <- (unique.bin.df$length / total.cov) * unique.bin.df[,seg.col.name]
    unique.bin.df$weighted.copy.ratio <- (unique.bin.df$length / total.cov) * unique.bin.df[,'copy.ratio']
    
    composite.sum <- sum(unique.bin.df$weighted.seg)
    composite.copy.ratio.sum <- sum(unique.bin.df$weighted.copy.ratio)
    composite.bin.final.df <- rbind(composite.bin.final.df, 
                                    data.frame(start.bin, end.bin, composite.sum, composite.copy.ratio.sum))
  }
  colnames(composite.bin.final.df) <- c("Start.bin", "End.bin", seg.col.name, 'copy.ratio')
  
  return(composite.bin.final.df)
}

# Function: generateGenomeBins
# Purpose:  Uses the UCSC ChromInfo file to generate bins of a specific size for later on
#   CN evaluation
# Input:  - chr.df <- ucsc.chrom.info dataframe containing Chr and Size
#         - bin.s <- The size of the bins you want to create
# Returns:  genome.binned: data-frame containing the genome binned into the set sizes
generateGenomeBins <- function(chr.df, bin.s){
  # Formats the ucsc chrom info dataframe in an easier to use way.
  chr.df <- chr.df[ chr.df$chrom %in% global.chrom,]
  rownames(chr.df) <- chr.df$chrom
  chr.df <- chr.df[,-1]
  
  # Creates a List containing binned intervals for ezch chromsome
  chr.bin.list <- list()
  for(each.chr in global.chrom){
    bin.start.pos <- 1
    bin.max.matrix <- matrix(ncol=2, nrow=0)
    while(bin.start.pos < max(chr.df[as.character(each.chr),'size'])){
      bin.end.pos <- bin.start.pos + bin.s
      if(bin.end.pos > max(chr.df[as.character(each.chr),'size'])){
        bin.end.pos <- max(chr.df[as.character(each.chr),'size'])
      }
      bin.max.matrix <- rbind(bin.max.matrix, c(bin.start.pos, bin.end.pos))
      bin.start.pos <- bin.end.pos + 1
    }
    chr.bin.list[[as.character(each.chr)]] <- bin.max.matrix
  }
  
  return(chr.bin.list)
}

# Function: addChrRect
# Purpose:  Adds blue transparent rectangles to indicate chromosome segments
# Input:  - type <- Default is a point plot, alternative is 'barplot' to adjust for segment size
# Returns:  NA
addChrRect <- function(type="point", max.ylim.sd=max.ylim.sd){
  box.start.num <- 0
  bp.scale.factor <- 1
  if(type %in% "barplot"){
    box.start.num <- 0.7  # Midpoint of the first box in a boxplot
    bp.scale.factor <- 1.2  # Each midpoint is 1.2 units aparts
  }
  
  # Adds blue labelled rectangles indicating chromosomal bins
  chrom.alt.count <- 1
  start.gen <- box.start.num
  # Takes the global chrom.bins.list variable for bins of chromosomes according to BP segment size
  for(each.chrom in 1:length(chrom.bins.list)){
    end.gen <- start.gen + (dim(chrom.bins.list[[each.chrom]])[1] * bp.scale.factor)
    if(chrom.alt.count %% 2 == 1){
      # Adds red rectangle
      rect(start.gen, -(max.ylim.sd), end.gen, max.ylim.sd, col=alpha("blue", 0.05), border=0, lwd=0)   
      
      # Adds text indicating Chr 
      text.pos <- ((end.gen - start.gen)/2) + start.gen
      text(text.pos, 0.05, labels=paste("Chr", chrom.alt.count, sep=""), cex=0.4)
    }
    start.gen <- end.gen
    
    chrom.alt.count <- chrom.alt.count + 1
  }
}


# Function: matchClAnno
# Purpose:  Given a list of CEL files, and a given matrix, it will identify which CEL files are annotated as the same cell line
# Input:  - x: A single row from the Annotation File
#         - match_val: either "match" or "mismatch" to identify what to look for and return
#         - total.matrix: The X-by-X matrix containing all the CELs to look up
# Returns:  A List containing either the matching values from the matrix or the mismatching values
matchClAnno <- function(x, match_val, total.matrix){
  match.values <- c()
  
  cl.fn <- as.vector(x[fn.headers])
  cl.fn <- cl.fn[!is.na(cl.fn)]
  rm.fn <- ""
  
  # Checks to make sure all the reported cell.line filenames from fn.headers are in cl.fn
  # removes ones that are not in that list.
  if(length(cl.fn[!(cl.fn %in% colnames(total.matrix))]) > 0){
    rm.fn <- cl.fn[!(cl.fn %in% colnames(total.matrix))]
    #print(paste("Not found: ", cl.fn[!(cl.fn %in% colnames(total.matrix))],
    #            "  --  ", "unique cell id: ", x['unique.cellid'], sep=""))
  } 
  
  # Checks to see if the number of non-matching cl.fn matches the length of cl.fn
  if(length(cl.fn[!(cl.fn %in% colnames(total.matrix))]) == length(cl.fn)){
    #print(paste("None of the dataset filenames were found in the match matrix:\n ", x['unique.cellid'], sep=""))
    match.values <- -1
    cl.fn <- 'NA'
    
    # Finds all the matching or mismatching based on the annotated matching cl.fn list
  } else {
    if(match_val == "match"){
      cl.fn <- cl.fn[cl.fn %in% colnames(total.matrix)]
      match.values <- total.matrix[cl.fn, cl.fn]
      if((class(match.values) == "numeric") && (length(match.values) == 1)){
        match.values <- as.matrix(match.values)
      }
      
      # Sets colnames and rownames to dataset names & cell line name for matching annotations
      cl.fn <- as.vector(x[fn.headers])
      # Removes na's and files not in the total.matrix
      fn.hd <- gsub(".filename", "", fn.headers[!(cl.fn %in% rm.fn) & !is.na(cl.fn)])
      fn.hd <- paste(fn.hd, x['unique.cellid'], sep="_")
      colnames(match.values) <- fn.hd
      rownames(match.values) <- fn.hd
      
    } else if(match_val == "mismatch"){
      cl.fn <- cl.fn[(cl.fn %in% colnames(total.matrix))]
      cl.mismatch.fn <- colnames(total.matrix)[!(colnames(total.matrix) %in% cl.fn)]
      if(!is.na(table(total.matrix[cl.fn,cl.mismatch.fn] %in% NA)["TRUE"])){
        print(cl.fn)
      }
      match.values <- total.matrix[cl.fn, cl.mismatch.fn]
    } 
  }
  
  if(class(match.values) == "numeric"){
    match.values <- t(match.values)
    rownames(match.values) <- cl.fn
  }
  
  return(as.matrix(match.values))
  
}

# Function: generateEasyClNames
# Purpose:  Takes a CEL name and finds its corresponding row and column number in the annotation file.  Will then reformat in an easy to read identifier
# Input:  - each.cl.name: CEL name
# Returns:  a String containing the renamed CEL file
generateEasyClNames <- function(each.cl.name){
  new.cl.name <- c()
  # Identifies the Row and Column in the cell-line annotation file that contains the variable
  row.index <- which(apply(cell.line.anno[,fn.headers] == each.cl.name, 1, any))
  col.index <- which(apply(cell.line.anno[,] == each.cl.name, 2, any))
  
  # If the cell line name is found in the annotation file, format it
  if(length(row.index) > 0  & length(col.index) > 0){
    cl.db <- sub(".filename.*", "", colnames(cell.line.anno[col.index]))
    if(length(grep("\\.y$", colnames(cell.line.anno[col.index]))) > 0){
      cl.db <- paste(cl.db, ".y", sep="")
    } else if (length(grep("\\.x$", colnames(cell.line.anno[col.index]))) > 0){
      cl.db <- paste(cl.db, ".x", sep="")
    }
    unique.cl <- cell.line.anno[row.index, 'unique.cellid']
    new.cl.name <- paste(cl.db, "_", unique.cl, sep="")

    # If cell line is not found, use the given name
  }  else {
    write(paste("Error: Trouble finding the name ", each.cl.name, " in annotation file!", sep=""), stderr())
    new.cl.name <- each.cl.name
  }
  return(new.cl.name[1])
}

# Function: relabelMatrix
# Purpose:  Simply renames the columns and rows for a given matrix
# Input:  - input.matrix: X-by-X matrix to relabel
#         - labels: the labels in a vector format
# Returns:  relabeled matrix
relabelMatrix <- function(input.matrix, labels){
  colnames(input.matrix) <- labels
  rownames(input.matrix) <- labels
  return(input.matrix)
}

# Function: generateMatchHistograms
# Purpose:  Generates the Match/Mismatch Histograms for a given X-by-X matrix
# Input:  - input.matrix: X-by-X matrix to plot 
#         - x.range: range of the x axis
#         - y.range: range of the y axis
#         - input.title: main title of the plot
#         - input.xlabel: x-axis label of the plot
#         - out.file: name of the output file
# Returns:  NA
generateMatchHistograms <- function(input.matrix, x.range, y.range, input.title, input.xlabel, out.file, num.breaks){
  match.list <- apply(cell.line.anno, 1, function(x) matchClAnno(x, "match", input.matrix))
  names(match.list) <- as.vector(cell.line.anno$unique.cellid)
  mismatch.list <- apply(cell.line.anno, 1, function(x) matchClAnno(x, "mismatch", input.matrix))
  names(mismatch.list) <- as.vector(cell.line.anno$unique.cellid)
  
  p1 <- hist(unlist(match.list), breaks=num.breaks)
  p2 <- hist(unlist(mismatch.list), breaks=num.breaks)
  pdf(file=out.file)
  plot(p1, xlim=x.range, ylim=y.range, col=rgb(0,0,1,1/4), main=input.title, xlab=input.xlabel)
  plot(p2, xlim=x.range, ylim=y.range, col=rgb(1,0,0,1/4), add=T)
  dev.off()
}

# Function: getHelp
# Purpose:  How to use this script help guide
# Input:  None  
# Returns:  String Vector: Statement of how to use the script
getHelp <- function(){
  help <- paste("###############################",
                "#       Tool Name             #",
                "###############################",
                "Purpose: ",
                "Usage: ",
                "Arguments: ", 
                sep="\n")
  
  return(help.message);
}


#########################################################################################################
#           Main
#########################################################################################################

setwd(cn.dir)
chrom.bins.list <- generateGenomeBins(ucsc.chrom.df ,global.bin.size)
total.bins <- sum(sapply(chrom.bins.list, function(x) length(x)))
print(paste("Total bins: ", total.bins, sep=""))
print(paste("Max Diff: ", max.diff, sep=""))

load(cl.anno.file)
#Allows for easy generation of Filenames later on
cell.line.anno$gdsc.filename.x <- gsub(".cel", ".CEL", cell.line.anno$gdsc.filename.x)
cell.line.anno$gdsc.filename.y <- gsub(".cel", ".CEL", cell.line.anno$gdsc.filename.y)

cn.files <- list.files(cn.dir, pattern="*segtab.txt")

if(process.data == 1){  
  # All against All CNV identity comparison
  #Cycles through each column header (hscr.a1, hscr.a2, or seg_means)
  for(seg.header in header.list){
    print(paste(Sys.time(), ": Starting analysis on ", seg.header, "...", sep=""))
    setwd(cn.dir)
    rdata.path <- file.path(output.dir, seg.header, rdata.dir)
    dir.create(rdata.path, showWarnings=F, recursive=TRUE)
    
    ##################
    ###### Reads in each segtab for the corresponding seg.header and stores the data frame in a Rdata file to load later
    if(create.df == 1){
      #cn.files <- list.files(cn.dir, pattern="*segtab.txt")
      print(paste("Spreading workload across ", num.of.cores, "...", sep=""))
      cl <- makeCluster(num.of.cores, outfile="", type="SOCK")
      registerDoSNOW(cl)
      
      foreach(each.cn=cn.files) %dopar% {
        require(intervals)
        out.fn <- gsub(".segtab.txt", paste(".", seg.header, ".Rdata", sep=""), each.cn)
        if(!file.exists(file.path(rdata.path, out.fn))){
          print(paste("Reading in: ", each.cn, sep=""))
          cn.file.df <- read.csv(each.cn, header=TRUE, sep="\t", comment.char="#")
          # Populates the empty bins with CN segment data
          cn.binned.df <- populateGenomeBins(cn.file.df[,c(chrom.header, start.header, end.header, 'copy.ratio', seg.header)], 
                                             chrom.bins.list, 
                                             seg.header)
          # Check for NA values in seg.header and copy.ratio and replace with missing values
          cn.binned.df[which(is.na(cn.binned.df[,'copy.ratio'])),c(seg.header, 'copy.ratio')] <- c(0,0)
          if(median(cn.binned.df[,'copy.ratio']) > 0){
            cn.binned.df[,seg.header] <- (cn.binned.df[,seg.header] / median(cn.binned.df[,'copy.ratio']))
          } else {
            # WHAT THE HELL DO I DO IF THE MEDIAN IS 0!?!?!
            write(paste("Error: ", each.cn, " has a median of 0.  Normaliation aborted!", sep=""), stderr())
          }
          save(cn.binned.df, file=file.path(rdata.path, out.fn))
        } else {
          print(paste("Previous record found and overwrite is not enabled: ", out.fn, sep=""))
        }
      }
      stopCluster(cl)
      # or uses a predesignated directory to load the existing Rdata files
    } else {
      if(existing.rdata.dir != ""){
        rdata.path <- existing.rdata.dir
      }
    }
    # A list containing all the Rdata files
    rdata.files <- list.files(rdata.path, pattern=paste("*\\.", seg.header, "\\.Rdata", sep=""))
    print(paste(Sys.time(), ": ", length(rdata.files), " .Rdata files are being compared...", sep=""))
    
    
    
    
    
    
    ##################
    ###### Generates the Squared-Difference Values
    print(paste(Sys.time(), ": Generating comprehensive binned-CN list: ", sep=""))
    total.rdata.list <- list()
    for(each.rdata in rdata.files){
      load(file.path(rdata.path, each.rdata))
      total.rdata.list[[each.rdata]] <- cn.binned.df     #Commented out and regenerated later after clust.list runs
    }
    
    
    print(paste(Sys.time(), ": Starting SSD analysis: ", sep=""))
    if(!file.exists(file.path(output.dir, seg.header, "sd_calc", "ssd.matrix.R"))){
      ssd.ptm <- proc.time()
      dir.create(file.path(output.dir, seg.header, "sd_calc"), showWarnings=FALSE, recursive=TRUE)
      ssd.score.df <- data.frame(matrix(c(rep(1, length(total.rdata.list))), nrow=1,
                                        dimnames=list(c(), names(total.rdata.list))),
                                 check.names=FALSE, stringsAsFactors=FALSE)   # Creates blank ssd matrix
      for(rdata.f1 in names(total.rdata.list)){
        ptm <- proc.time()
        val <- grep(rdata.f1, names(total.rdata.list))
        print(paste("Comparing ", rdata.f1, sep=""))  #Cell line name and cn.binned.df stored in that list cluster
        subset.comp <- lapply(total.rdata.list[c(val:length(total.rdata.list))], 
                              function(x) getSsdScoreTally(total.rdata.list[[rdata.f1]], x, 
                                                           seg.header, total.bins, max.diff))
        matrix.ssd.row <- as.vector(unlist(subset.comp))
        names(matrix.ssd.row) <- names(total.rdata.list)[c(val:length(total.rdata.list))]
        ssd.score.df <- smartbind(ssd.score.df, matrix.ssd.row)
        
        save(ssd.score.df, file=file.path(output.dir, seg.header, "sd_calc", "ssd.matrix.R"))
        print(proc.time() - ptm)
      }
      ssd.score.df <- ssd.score.df[-1,]
      rownames(ssd.score.df) <- names(total.rdata.list)
      ssd.score.df <- fillSimilarityMatrix(ssd.score.df)  #Fills other half of triangle
      save(ssd.score.df, file=file.path(output.dir, seg.header, "sd_calc", "ssd.matrix.R"))
      print(paste(Sys.time(), ": Finished SSD analysis in: ", proc.time() - ptm,  sep=""))
    } else {
      print(paste(Sys.time(), ": SD-Matrix already found, loading previous iteration. ", sep=""))
      load(file.path(output.dir, seg.header, "sd_calc", "ssd.matrix.R"))
    }
    
    # Plots SD plots for all files in the plot.list (args[8])
    plots.path <- file.path(output.dir, seg.header, plots.dir)
    if(length(plot.list) > 0){
      plotSd(file.l=plot.list, comp.all=FALSE, total.rdata.list=total.rdata.list,
             seg.header=seg.header, plots.path=plots.path)
    }
    
    
    
    
    ##################
    ###### Generates the Squared-Difference Values
    print(paste(Sys.time(), ": Finished generating bulk data and SD plots.", sep=""))
    print(paste(Sys.time(), ": Starting post-processing and outputting of RData frames", sep=""))
    ### Post-Processing Data Formatting and Saving
    dir.create(file.path(output.dir, seg.header, "data_matrix"), showWarnings=FALSE, recursive=TRUE)
    data.matrix.dir <- file.path(output.dir, seg.header, "data_matrix")
    save(cn.files, file=file.path(data.matrix.dir, out.cn.files))
    save(total.rdata.list, file=file.path(data.matrix.dir, out.total.rdata))
    
    # Relabelling Matrices as .CEL files
    total.cn.matrix <- as.matrix(ssd.score.df)
    total.cn.matrix <- relabelMatrix(total.cn.matrix, gsub(".segtab.txt", ".CEL", cn.files))
    rm(ssd.score.df)
    
    # Creates the Spearman R Correlation Matrix
    total.matrix <- do.call("cbind", lapply(total.rdata.list, function(x) x[,seg.header]))
    spear.corr.matrix <- cor(total.matrix, method = "spearman")
    spear.corr.matrix <- relabelMatrix(spear.corr.matrix, gsub(".segtab.txt", ".CEL", cn.files))
    spear.corr.simple.matrix <- relabelMatrix(spear.corr.matrix, 
                                              unlist(lapply(colnames(spear.corr.matrix), function(x) generateEasyClNames(x))))
    #Generates a eucldiean distance plot
    #eucl.dist.matrix <- dist(total.matrix)
    #eucl.dist.matrix <- relabelMatrix(eucl.dist.matrix, gsub(".segtab.txt", ".CEL", cn.files))
    
    #Rdata file outputs
    save(total.cn.matrix, file=file.path(data.matrix.dir, out.cn.sse))
    save(spear.corr.matrix, file=file.path(data.matrix.dir, out.spearman.mat))
    save(spear.corr.simple.matrix, file=file.path(data.matrix.dir, paste("simple", out.spearman.mat, sep=".")))
    #save(eucl.dist.matrix, file=out.eucl.mat)
    
    print(paste(Sys.time(), ": Finished processing the data for: ", seg.header, sep=""))
    print(paste(Sys.time(), ": Generating heatmap for: ", seg.header, sep=""))
    #### PLOT #5: #Plots the heatmap for the -log of normalized sum of cubed error for CN segment comparison
    # Plot 5a) No threshold applied to the data
    # Heatmap output
    print("Generating 'Plot 5': Heatmap of the -Log(Norm[SCD])...")
    if(capabilities()['png']){
      out.file <- paste(seg.header, ".heatmap_ssd.png", sep="")
      png(file=file.path(data.matrix.dir, out.file), width=20, height=20,units="in",res=800)
    } else {
      out.file <- paste(seg.header, ".heatmap_ssd.pdf", sep="")
      pdf(file=file.path(data.matrix.dir, out.file), width=20, height=20)
    }
    my_palette <- colorRampPalette(c("red","black", "green"))(n = 10000)
    heatmap.2(spear.corr.matrix, margin=c(10,10), trace="none",col=my_palette,
              dendrogram="none",
              lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ),
              key=FALSE,
              cexRow=0.05,cexCol=0.05)
    dev.off()
    print(paste(Sys.time(), ": Finished heatmap output.  Completed run of : ", seg.header, sep=""))
  }# End of seg.header 
}


# #####
# #TEST SET:
# #output.dir <- '/Users/rquevedo/Desktop/bhk_lab/results/cnv_fingerprinting/50k_subset'
# #########################################################################################################
# # Input List
# if(generate.output == 1 & ignore.process == 0){
#   complete.matrix.list <- list()
#   matrix.count <- 1
#   
#   for(seg.header in header.list){
#     # Load the seg file for the header
#     setwd(file.path(output.dir, seg.header, 'data_matrix'))
#     load(out.cn.sse)
#     load(out.thresh.cn.see)
#     load(out.cn.mat)
#     load(out.cn.files)
#     load(out.total.rdata)
#     load(out.spearman.mat)
#     load(out.tmat.hscr)
#     
#     #Package the matrices and vectors into a List
#     seg.header.list <- list()
#     seg.header.list[['cn.files']] <- cn.files
#     seg.header.list[['total.rdata.list']] <- total.rdata.list
#     seg.header.list[['total.cn.matrix']] <- total.cn.matrix
#     seg.header.list[['total.thresh.cn.matrix']] <- total.thresh.cn.matrix
#     seg.header.list[['total.cn.list.matrix']] <- total.cn.list.matrix
#     seg.header.list[['spear.corr.matrix']] <- spear.corr.matrix
#     seg.header.list[['total.matrix']] <- total.matrix
#     complete.matrix.list[[matrix.count]] <- seg.header.list
#     matrix.count <- matrix.count + 1
#     
#     # Create the Total-HSCR Data set from the A- and B- allele HSCR
#     if(matrix.count == 3){
#       total.cn.matrix <- complete.matrix.list[[1]][['total.cn.matrix']] + complete.matrix.list[[2]][['total.cn.matrix']]
#       total.thresh.cn.matrix <- complete.matrix.list[[1]][['total.thresh.cn.matrix']] + complete.matrix.list[[2]][['total.thresh.cn.matrix']]
#       total.cn.matrix <- complete.matrix.list[[1]][['total.cn.matrix']] + complete.matrix.list[[2]][['total.cn.matrix']]
#       total.matrix <- complete.matrix.list[[1]][['total.matrix']] + complete.matrix.list[[2]][['total.matrix']]
#       # Creates the Spearman R Correlation Matrix
#       spear.corr.matrix <- cor(total.matrix, method = "spearman")
#       spear.corr.matrix <- relabelMatrix(spear.corr.matrix, gsub(".segtab.txt", ".CEL", cn.files))
#       
#       seg.header.list <- list()
#       seg.header.list[['cn.files']] <- cn.files
#       seg.header.list[['total.rdata.list']] <- total.rdata.list
#       seg.header.list[['total.cn.matrix']] <- total.cn.matrix
#       seg.header.list[['total.thresh.cn.matrix']] <- total.thresh.cn.matrix
#       seg.header.list[['total.cn.list.matrix']] <- total.cn.list.matrix
#       seg.header.list[['spear.corr.matrix']] <- spear.corr.matrix
#       seg.header.list[['total.matrix']] <- total.matrix
#       complete.matrix.list[[matrix.count]] <- seg.header.list
#     }
#   }
#   
#   matrix.count <- 1
#   header.list[[3]] <- 'total.cn'
#   for(each.matrix in complete.matrix.list){
#     seg.header <- header.list[[matrix.count]]
#     # Load the Correct Files to Process
#     cn.files <- each.matrix[['cn.files']]
#     total.rdata.list  <- each.matrix[['total.rdata.list']]
#     total.cn.matrix  <- each.matrix[['total.cn.matrix']]
#     total.thresh.cn.matrix  <- each.matrix[['total.thresh.cn.matrix']]
#     total.cn.list.matrix  <- each.matrix[['total.cn.list.matrix']]
#     spear.corr.matrix  <- each.matrix[['spear.corr.matrix']]
#     total.matrix  <- each.matrix[['total.matrix']]
# 
#     
#     #########################################################################################################
#     #         Output
#     ######################################################################################################### 
#     dir.create(file.path(output.dir, "total_analysis"), showWarnings=FALSE, recursive=TRUE)
#     setwd(file.path(output.dir, "total_analysis"))
#     
#     #### PLOT #4: #Histogram of matched and mismatched CEL files and their SCD values
#     # Plot 4a) No threshold applied to the data
#     print("Generating 'Plot 4': Histograms showing the distribution of match vs mismatched CEL CN-SCD values...")
#     generateMatchHistograms(total.cn.matrix, c(0,4), c(0,1000), 
#                             "-Log Sum of Cubed Error of Binned Copynumber", 
#                             "-Log CN-binned SCD", 
#                             paste(seg.header, ".cn.match_hist.pdf", sep=""),
#                             1000)
#     
#     # Plot 4b) Local minima threshold between an acceptable range applied to the data
#     generateMatchHistograms(total.thresh.cn.matrix, c(0,10), c(0,1000), 
#                             "-Log Sum of Cubed Error of Binned Copynumber - Threshold", 
#                             "-Log CN-binned SCD", 
#                             paste(seg.header, ".thresh.cn.match_hist.pdf", sep=""),
#                             1000)
#     
#     # Plot 4c) Values based on a spearman correlation matrix
#     generateMatchHistograms(spear.corr.matrix, c(-1,1), c(0,1000), 
#                             "Spearman's R Histogram", 
#                             "Spearman's R", 
#                             paste(seg.header, ".spearman.cn.match_hist.pdf", sep=""),
#                             200)
#     
#     # Plot 4d) Values based on a composite spearman R and Norm-SCD score
#     comp.matrix <- (spear.corr.matrix * 0.5) + ((1-exp(-total.thresh.cn.matrix)) * 0.5)
#     comp.matrix <- relabelMatrix(comp.matrix,  gsub(".segtab.txt", ".CEL", cn.files))
#     
#     generateMatchHistograms(comp.matrix, c(min(comp.matrix),0), c(0,10), 
#                             "Composite SCD & Spearman's R Score", 
#                             "Composite Score", 
#                             paste(seg.header, ".composite.cn.match_hist.pdf", sep=""),
#                             200)
#     
#     
#     # Convert the names to easier to read format
#     raw.cel.names <- gsub(".segtab.txt", ".CEL", cn.files)
#     easy.cel.names <- unlist(lapply(colnames(matrix.perc.match), function(x) generateEasyClNames(x)))
#     
#     total.cn.matrix <- relabelMatrix(total.cn.matrix, easy.cel.names)
#     total.thresh.cn.matrix <- relabelMatrix(total.thresh.cn.matrix, easy.cel.names)
#     matrix.perc.match <- relabelMatrix(matrix.perc.match, easy.cel.names)
#     
#     
#     #### PLOT #5: #Plots the heatmap for the -log of normalized sum of cubed error for CN segment comparison
#     # Plot 5a) No threshold applied to the data
#     # Heatmap output
#     print("Generating 'Plot 5': Heatmap of the -Log(Norm[SCD])...")
#     out.file <- paste(seg.header, ".heatmap_ssd.png", sep="")
#     png(file=out.file, width=20, height=20,units="in",res=800)
#     my_palette <- colorRampPalette(c("red","black", "green"))(n = 10000)
#     heatmap.2(spear.corr.matrix, margin=c(10,10), trace="none",col=my_palette,
#               dendrogram="none",
#               lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ),
#               key=FALSE,
#               cexRow=0.05,cexCol=0.05)
#     dev.off()
#     
#     # Plot 5b) Local minima threshold between an acceptable range applied to the data
#     # Threshold Heatmap output
#     out.file <- paste(seg.header, ".thresh.heatmap_ssd.png", sep="")
#     png(file=out.file, width=20, height=20,units="in",res=800)
#     my_palette <- colorRampPalette(c("red","white", "blue"))(n = 10000)
#     heatmap.2(matrix.perc.match, margin=c(10,10), trace="none",col=my_palette,
#               dendrogram="none",
#               lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ),
#               key=TRUE,
#               cexRow=0.1,cexCol=0.1)
#     dev.off()
#     
#     #### PLOT #6: #Plots the heatmap for the -log of normalized sum of cubed error for CN segment comparison
#     print("Generating 'Plot 6': Correlogram using Spearmans correlation between each cell lines CN profile...")
#     # Plot 6a) Spearman Correlogram using Shade and Ellipse display
#     out.file <- paste(seg.header, ".spearman_corrgram.ellipse.png", sep="")
#     png(file=out.file, width=20, height=20,units="in",res=800)
#     corrgram(total.matrix, cor.method="spearman",
#              order=TRUE, 
#              lower.panel=panel.shade, 
#              upper.panel=panel.ellipse, 
#              text.panel=panel.txt,) 
#     dev.off()
#     
#     print(paste("Successfully outputted all plots for ", seg.header, sep=""))
#     matrix.count <- matrix.count + 1
#   }
# }
#   
#   
