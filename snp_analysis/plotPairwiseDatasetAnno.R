.libPaths(c(.libPaths(), "/mnt/work1/users/pughlab/bin/r_lib/library/3.1"))
###############################################################
#
#                 SNP Dataset Heatmap And Visualization
#
# Author: Rene Quevedo
# Date Created: Dec-07-2015
###############################################################
# Function:
#   - Creates a heatmap of percent identity.
#   - Creates a density plot showing mismatched cell lines based
#         on annotations
###############################################################
library(gplots)
library(gtools)
library(scales)
library(reshape)
library(ggplot2)
source('/Users/rquevedo/git/cnv_fingerprint/v2.0/func/processAnnotations.R')
load('/Users/rquevedo/git/cnv_fingerprint/v2.0/data/merged_annotations.Rdata')

data.type <- 'snp' # either cnv or snp
low.match.threshold <- c("cnv"=0.6,
                         "snp"=0.8)

mismatch.col <- 'magenta'
match.col <- 'brown'

snp.matrix.dir <- '~/Onedrive/bhk_lab/data/snp_fingerprinting'
snp.matrix.dir <- '/data/snp/'
snp.matrix.file <- 'complete_snp_matrix.allSnps.R'

global.curation.dir <- '/data/ref/'
out.rdata <- 'merged_annotations.Rdata'

global.wd <- '/results/snp/'
dir.create(global.wd, recursive = TRUE, showWarnings = FALSE)


coi <- NA
coi <- c("NCI-H1869", "HuP-T4", 'Ishikawa (Heraklio) 02 ER-') # Cell lines of interest instead of concordance
coi <- c("NCI-H1869", "HuP-T4", 'NCI-H23') # Cell lines of interest instead of concordance

coi.col <- 'red'

load(file.path(global.curation.dir, out.rdata))
load(file.path(snp.matrix.dir, snp.matrix.file))
matrix.perc.match <- t.matrix.perc.match

###############################
#         Functions
###############################
plotNoncordantCl <- function(x, title, target, binarize=FALSE,
                             nonmatch.threshold=0.7, ambigious.threshold=0.8){
  myPalette <- colorRampPalette(c(mismatch.col, "white"), bias=1, interpolate="linear")

  x.m <- melt(as.matrix(x))
  colnames(x.m) <- c("cell.line", "dataset", "concordance")
  x.m$dataset <- factor(x.m$dataset, levels=target)

  if(binarize) {
    x.m$bin <- cut(x.m$concordance,
                   breaks=c(0, nonmatch.threshold,
                            ambiguous.threshold, 1))
    sc <- scale_fill_manual(values = myPalette(3), na.value="grey50")
  } else {
    x.m$bin <- round(x.m$concordance,2)
    sc <- scale_fill_gradientn(colours = myPalette(100), na.value="grey50")
  }

  p <- ggplot(x.m, aes(dataset, cell.line)) +
    geom_tile(aes(fill = bin), colour = "black") +
    sc

  p + theme_bw() +
    ggtitle(title) +
    labs(x="", y="") +
    theme(legend.position = "right",
          axis.ticks = element_blank(),
          axis.ticks.length=unit(0.85, "cm"),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 0, colour = "black")) +
    coord_fixed(ratio=1)
}



# Function: addAnnoLines
# Purpose:  Adds annotations and lines to given spoints
# Input:   anno.points <- a named vector of percent identities
# Returns:
addAnnoLines <- function(anno.points, max.height=0.75, x.pos=1, dist.y=0.05, dir='right'){
  anno.points <- low.match
  anno.points <- anno.points[order(anno.points, decreasing = TRUE)]
  anno.names <- names(anno.points)
  if(dir == 'right'){
    x.pos.1 <- (x.pos + 1)
    x.pos.2 <- (x.pos + 2)
    x.text <- (x.pos + 2.1)
  } else if (dir == 'left'){
    x.pos.1 <- (x.pos - 1)
    x.pos.2 <- (x.pos - 2)
    x.text <- (x.pos - 2.1)
  } else {
    stop("Requires a 'right' or 'left' direction")
  }

  # Diagonal line from points
  segments(x0=rep(x.pos, length(anno.points)),
           y0=as.numeric(anno.points),
           x1=rep(x.pos.1, length(anno.points)),
           y1=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y))
  # Horizontal line from diagonal to name
  segments(x0=rep(x.pos.1, length(anno.points)),
           y0=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y),
           x1=rep(x.pos.2, length(anno.points)),
           y1=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y))
  # Name of Segments
  text(x=rep(x.text, length(anno.points)),
       y=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y),
       labels=anno.names,
       adj=if(dir == 'right') 0 else 1)
}





###############################
#           Main
###############################

setwd(global.wd)
#Load and reformat/rename the files
rownames(matrix.perc.match) <- gsub(".cel$", ".CEL", rownames(matrix.perc.match))
colnames(matrix.perc.match) <- gsub(".cel$", ".CEL", colnames(matrix.perc.match))
cell.line.anno$gdsc.filename.x <- gsub(".cel$", ".CEL", cell.line.anno$gdsc.filename.x)
cell.line.anno$gdsc2.filename.y <- gsub(".cel$", ".CEL", cell.line.anno$gdsc2.filename.y)
col.index <- which(colnames(cell.line.anno) %in% c('pfizer2.filename.y', 'pfizer3.filename', 'gdsc2.filename.y'))
colnames(cell.line.anno)[col.index] <- c("pfizer2.filename.y", "pfizer3.filename", "gdsc2.filename.y")

# Match CEL to cell lines and generate match/mismatches
low.match.list <- list()
ds.type <- c('cgp.filename', 'ccle.filename', 'gdsc.filename.x', 'pfizer.filename.x',
             'pfizer.filename.y', 'pfizer.filename', 'gdsc2.filename.y')
for(ds1 in ds.type){
  low.match.df <- data.frame(matrix(ncol=1, nrow=0))
  low.match.list[[ds1]] <- low.match.df
  for(ds2 in ds.type){
    print(ds2)
    if(!ds1 == ds2){
      rm(match.list, mismatch.list)

      fn.headers <- c(ds1, ds2)  # colnames in cell.line.anno for CGP and CCLE filenames
      match.list <- apply(cell.line.anno, 1, function(x) matchClAnno(x, "match", matrix.perc.match, fn.headers, verbose=FALSE))
      names(match.list) <- as.vector(cell.line.anno$unique.cellid)
      mismatch.list <- apply(cell.line.anno, 1, function(x) matchClAnno(x, "mismatch", matrix.perc.match, fn.headers))
      names(mismatch.list) <- as.vector(cell.line.anno$unique.cellid)


      ###############################
      #           Visualization
      ###############################

      # Creates a heatmap of all against all concordance scores
      # print("Creating heatmap...")
      #
      # sunset.col <- colorRampPalette(c("black","red", "orange", "yellow", "white"))(100)
      # print("Creating a heatmap of the match.")
      # png(file=out.file, width=20, height=20,units="in",res=500)
      # heatmap.2(matrix.perc.match, margin=c(50,50), trace="none",col=sunset.col,
      #           lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ),
      #           cexRow=0.12,cexCol=0.12)
      # dev.off()
      #
      # print("Completed SNP Fingerprinting")
      # date()

      #Create match for each cell line density and scatterplot
      match.v <- unlist(match.list, use.names=TRUE)
      match.v <- match.v[which(match.v != 1)]   # cell lines that match to themselves, diagonal
      match.v <- match.v[which(match.v != -1)]  #
      if(length(match.v) > 0){
        match.v <- match.v[-seq(1, length(match.v), by=2)] #Remove duplicated entries
        names(match.v) <- gsub("[1234]$", "", names(match.v), perl=TRUE) #Remove index number of matrix
        dens.mm <- density(unlist(mismatch.list))
        dens.match <- tryCatch(density(match.v),
                               error = function(e) {density(rep(match.v, 2))})
        #low.match.threshold[data.type] <- quantile(match.v, c(0.05))
        #Create mismatch for each cell line density and scatterplot
        #mismatch.v <- unlist(mismatch.list, use.names=TRUE)

        #mismatch.v <- unlist(mismatch.list, use.names=TRUE)
        quant.match.cutoff <- quantile(match.v, 0.05)

        #http://www.r-bloggers.com/absolute-deviation-around-the-median/
        quant.match.cutoff <- median(match.v) - (3*mad(match.v, center=median(match.v), constant=1.4826, na.rm=TRUE))
        #constant=(1/quantile(match.v, 0.75))
        if(quant.match.cutoff < low.match.threshold[data.type]) {cutoff.val <- low.match.threshold[data.type] } else { cutoff.val <- quant.match.cutoff}
        if(cutoff.val > 0.9 | cutoff.val < 0.6){cutoff.val <- 0.9}



        #pdf(paste(paste(gsub(".filename(.+)?", "", fn.headers, perl=TRUE), sep="", collapse="_"), ".annoplot.pdf", sep=""))
        layout(matrix(c(rep(c(0,0,2,2,1,0), 9),
                        c(0,0,0,3,3,0)), ncol=6, byrow=TRUE))
        par(mar=c(0, 0, 4.1, 2.1))
        #par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
        plot(dens.match$y, dens.match$x, type="l", col=match.col,
             ylim=c(if(data.type=='cnv') -1 else 0,1),
             xlab="", ylab="",
             yaxt='n', xaxt='n')
        lines((dens.mm$y * (max(dens.match$y) / max(dens.mm$y))),
              dens.mm$x,, col=mismatch.col)
        # Add in coloured shades
        with(dens.match, polygon(x=dens.match$y,
                                 y= dens.match$x,
                                 col=alpha(match.col, 0.5)))
        with(dens.mm, polygon(x=(dens.mm$y * (max(dens.match$y) / max(dens.mm$y))),
                              y= dens.mm$x,
                              col=alpha(mismatch.col, 0.5)))

        abline(h = low.match.threshold[data.type], col="grey", lty=5)
        abline(h = cutoff.val, col="dimgrey", lty=3)

        par(mar=c(0, 4.1, 4.1, 0))
        plot(rep(1, length(unique(match.v))), unique(match.v), col=alpha(match.col, 0.1),
             ylim=c(if(data.type=='cnv') -1 else 0,1), xlim=c(0,10), las=2,
             xaxt='n', xlab='',
             ylab=if(data.type=='cnv') 'Spearman Correlation' else 'Concordance',
             main=paste(gsub(".filename(.+)?", "", fn.headers, perl=TRUE), sep="", collapse=" - "))

        text(x=3.1,
             y=1.0, labels = paste("n = ", length(match.v), sep=""),
             adj=0)

        # Add in high-match concordance "mismatch" points
        #high.mismatch <- mismatch.v[which(mismatch.v >= low.match.threshold[data.type])]
        #points(rep(1, length(high.mismatch)), high.mismatch, col=mismatch.col)

        # Add in low-match concordance "matching" points
        low.match <- match.v[which(match.v < cutoff.val)]
        if(!all(is.na(coi))) low.match <- match.v[which(names(match.v) %in% coi)]

        points(rep(1, length(low.match)),  low.match,
               pch=if(!all(is.na(coi))) 16 else 1,
               col=if(!all(is.na(coi))) coi.col  else match.col)
        if(length(low.match) > 0){
          addAnnoLines(low.match,
                       max.height=(cutoff.val - 0.05),
                       dist.y=0.03,
                       x.pos=1,
                       dir='right')
          abline(h = low.match.threshold[data.type], col="grey", lty=5)
          abline(h = cutoff.val, col="dimgrey", lty=3)
        }

        # Plot the legend
        par(mar=c(0,0,0,0))
        plot(0, type='n', xlim=c(0,10), ylim=c(0,10), axes=FALSE, ylab='', xlab='')
        legend(0.6,
               10, # places a legend at the appropriate place
               c('matching cell lines','mismatching cell lines'), # puts text in the legend
               lty=c(1,1), # gives the legend appropriate symbols (lines)
               lwd=c(2.5,2.5),col=c('brown',mismatch.col)) # gives the legend lines the correct color and width
        dev.off()


        low.match.list[[ds1]] <- t(smartbind(t(low.match.list[[ds1]]), t(data.frame(low.match))))
      } else {
        # If there is no matches found between the two datasets
        low.match <- NA
        low.match.list[[ds1]] <- t(smartbind(t(low.match.list[[ds1]]), t(data.frame(low.match))))
      }

    } else {
      # if ds2 = ds1
      low.match <- NA
      low.match.list[[ds1]] <- t(smartbind(t(low.match.list[[ds1]]), t(data.frame(low.match))))
    }
  }
  low.match.list[[ds1]] <- low.match.list[[ds1]][,-1, drop=FALSE]
  colnames(low.match.list[[ds1]]) <- ds.type
}

pre.low.match.list <- low.match.list

## Find overlap between "low.match" cell lines and existence in other datasets, fills in matches
for(each.ds.1 in names(low.match.list)){
  ds1.index <- grep(each.ds.1, ds.type)[1]
  all.ds.2 <- colnames(low.match.list[[each.ds.1]])
  for(each.ds.2 in all.ds.2){
    ds2.index <- grep(paste("^", each.ds.2, "$", sep=""), ds.type)[1]

    # Gets "matching" annotated CEL files for ds.1 and ds.2
    cl.row.index <- match(rownames(low.match.list[[each.ds.1]]), cell.line.anno$unique.cellid)
    ds1.CEL.files <- as.character(cell.line.anno[cl.row.index, which(colnames(cell.line.anno) %in% ds.type[ds1.index])])
    ds2.CEL.files <- as.character(cell.line.anno[cl.row.index, which(colnames(cell.line.anno) %in% ds.type[ds2.index])])

    #Indexes based on the similarity matrix
    matrix.ds1.index <- match(ds1.CEL.files, rownames(matrix.perc.match))
    matrix.ds2.index <- match(ds2.CEL.files, colnames(matrix.perc.match))

    #Obtain a pairwise concordance score if it's there, else NA
    match.val <- c()
    for(each.val in 1:length(matrix.ds1.index)){
      match.val <- c(match.val,
                     matrix.perc.match[matrix.ds1.index[each.val], matrix.ds2.index[each.val]])
    }

    #Replace only the NA values that didnt have a low-match score
    match.val[match.val > 0.8] <- 1
    replace.index <- which(is.na(low.match.list[[each.ds.1]][,each.ds.2]))
    low.match.list[[each.ds.1]][replace.index, each.ds.2] <- match.val[replace.index]
  }
  colnames(low.match.list[[each.ds.1]]) <- gsub(".filename", "", colnames(low.match.list[[each.ds.1]]))
  colnames(low.match.list[[each.ds.1]]) <- gsub("pfizer$", "pfizer3", colnames(low.match.list[[each.ds.1]]))
  colnames(low.match.list[[each.ds.1]]) <- gsub("pfizer.y", "pfizer2", colnames(low.match.list[[each.ds.1]]))
  colnames(low.match.list[[each.ds.1]]) <- gsub("\\.[xy]", "", colnames(low.match.list[[each.ds.1]]))

}
names(low.match.list) <- gsub(".filename", "", names(low.match.list))
names(low.match.list) <- gsub("pfizer$", "pfizer3", names(low.match.list))
names(low.match.list) <- gsub("pfizer.y", "pfizer2", names(low.match.list))
names(low.match.list) <- gsub("\\.[xy]", "", names(low.match.list))
lapply(low.match.list, function(x) {colnames(x) <- names(low.match.list); return(x)})

pdf("non-concordant-ccl.ccle.all.pdf")
lapply(names(low.match.list), function(x) plotNoncorcodantCl(low.match.list[[x]], x, colnames(low.match.list[[x]])))
dev.off()



# Summary statistics for Cell lines with same annotation but at least 1 discordant pairwise match
x.df <- do.call("rbind", low.match.list)
y.df <- x.df[order(rownames(x.df)),]
rownames(y.df)[which(!duplicated(rownames(y.df)))]  # Num of unique cell lines affected

#Save environment
save(y.df, low.match.list, cell.line.anno, matrix.perc.match,
     file="/Users/rquevedo/git/snp_fingerprint/environments/plotPairwiseSnp.v2.Rdata")

# DELETE FROM HERE DOWN: Used to get concordance and correlation values for cell lines
x <- unlist(cell.line.anno[grep("^HUH-6-clone5$", cell.line.anno$unique.cellid), c("cgp.filename", "gdsc.filename.x")])
x <- gsub("cel", "CEL", x)
x <- x[!is.na(x)]

hscr.a1.mat[x,x]
hscr.a2.mat[x,x]
t.matrix.perc.match[x,x]
