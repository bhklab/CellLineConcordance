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
library(scales)
library(gtools)
library(reshape)
library(ggplot2)

data.type <- 'cnv' # either cnv or snp
low.match.threshold <- c("cnv"=0.6,
                         "snp"=0.8)
max.num.low.match <- 20
report.summary <- TRUE

mismatch.col <- 'magenta'
match.col <- 'brown'
mismatch.col <- "gray67"
match.col <- "gray9"
coi <- c("NCI-H1869", "HuP-T4", 'Ishikawa (Heraklio) 02 ER-') # Cell lines of interest instead of concordance
coi <- c("NCI-H1869", "HuP-T4", 'NCI-H23') # Cell lines of interest instead of concordance
coi.col <- "red"  # Colour for points matched to cell lines of interest

pdir <- '~/git/CellLineConcordance/cnv_analysis/'
global.curation.dir <- file.path(pdir, "data");
out.rdata <- 'merged_annotations.Rdata';

cnv.matrix.dir <- file.path(pdir, "data");
cnv.matrix.file <- file.path(cnv.matrix.dir, "hscra1.weightedCorrDf.Rdata")
load(cnv.matrix.file)
hscr.a1.mat <- weighted.corr.df
cnv.matrix.file <- file.path(cnv.matrix.dir, "hscra2.weightedCorrDf.Rdata")
load(cnv.matrix.file)
hscr.a2.mat <- weighted.corr.df
global.wd <- "~/git/CellLineConcordance/cnv_analysis/results/cnv"
dir.create(global.wd, recursive=TRUE)
setwd(global.wd)

load(file.path(global.curation.dir, out.rdata))

###############################
#         Functions
###############################
{
rmNaRows <- function(x){
  na.rows <- which(apply(x, 1, function(y) all(is.na(y))))
  x <- x[-na.rows,]
  return(x)
}

# title <- names(low.match.list)[1]
# x <- low.match.list[[title]]
# out.dir.nc <- file.path(global.wd, "allele_nc")

plotNoncorcodantCl <- function(x, title, out.dir.nc){
  x.a <- forceMatrix(rmNaRows(x[,grep("_a", colnames(x))]), x)
  x.b <- forceMatrix(rmNaRows(x[,grep("_b", colnames(x))]), x)

  if(dim(x.a)[1] == 0 | dim(x.b)[1] == 0 ){
    warning(paste("No discordance found for ", title, sep=""))
    return(NULL)
  }

  x.a.b.list <- list("a"=x.a, "b"=x.b)

  # Melt and a_ and b_ allele discordances
  for(allele in c("a", "b")){
    x.temp <- x.a.b.list[[allele]]
    #x.temp[is.na(x.temp)] <- 1
    x.temp <- (1-x.temp)
    # if(is.matrix(x.temp))
    colnames(x.temp) <- gsub("_[ab]$", "", colnames(x.temp))
    x.m <- melt(x.temp)
    colnames(x.m) <- c("cell.line", "dataset", "concordance")
    target <- c("ccle", "cgp", "gdsc", "gdsc2", "pfizer", "pfizer2", "pfizer3")
    x.m$dataset <- factor(x.m$dataset, levels=target)
    x.m$allele <- rep(allele, dim(x.temp)[1])
    x.a.b.list[[allele]] <- x.m
  }

  x.a.b.list <- matchMeltedAlleles(x.a.b.list)
  x.ab <- do.call("rbind", x.a.b.list)

  #Replace the pseudo-filled in 1's with the actual concordances after checking if they're in the same order
  real.val.v <- c(melt(x.a)$value, melt(x.b)$value)
  # real.val.v[is.na(real.val.v)] <- 1
  real.val.v <- (1-real.val.v)
  real.val.chk <- real.val.v
  real.val.chk[real.val.chk > 0] <- 1
  if(all(real.val.chk == x.ab$concordance, na.rm=TRUE)){
    x.ab$concordance  <- real.val.v
  }

  split.melt <- split(x.ab, x.ab[,'match'])
  for(each.match in names(split.melt)){
    # Sets flags to print a header (first page) or a legend (last page)
    if(each.match == names(split.melt)[1]) print.head <- 1 else print.head <- 0
    if(each.match == names(split.melt)[length(names(split.melt))]) print.legend <- 1 else print.legend <- 0

    # If only one allele df is present, generates the opposite allele with all NAs
    if(each.match %in% 'A-absent'){
      split.melt[[each.match]] <- rbind(genMissingAlleleNa(split.melt[[each.match]], "a"), split.melt[[each.match]])
    }
    if(each.match %in% 'B-absent'){
      split.melt[[each.match]] <- rbind(split.melt[[each.match]], genMissingAlleleNa(split.melt[[each.match]], "b"))
    }

    # Prints the subheatmaps
    pdf(file.path(out.dir.nc, paste(title, ".", each.match, ".nc-ccl.pdf", sep="")))
    print(paste("Outputing summary match: ", paste(title, ".", each.match, ".nc-ccl.pdf", sep="")))
    p <- ggplot(split.melt[[each.match]], aes(dataset, cell.line)) +
      geom_tile(aes(fill = concordance), colour = "black") +
      coord_fixed(ratio = 1) +
      scale_fill_gradient(low = "white",high = mismatch.col,  na.value="grey50", limits=c(0, 1.5))
    print(p + theme_bw() +
            ggtitle(if(print.head == 1) title else "") +
            labs(x="", y="") +
            theme(legend.position = if(print.legend == 1) "bottom" else "none",
                  strip.text.x=if(print.head == 0) element_blank() else element_text(),
                  axis.ticks = element_blank(),
                  axis.text.x = element_text(size = 10, angle = 90, hjust = 0, colour = "black")) +
            facet_grid(match ~ allele))
    dev.off()
  }
}

forceMatrix <- function(x.temp, ref.x){
  if(!is.matrix(x.temp)){
    match.l <- apply(ref.x[,names(x.temp)], 1, function(y) all(x.temp %in% y))
    cl.name <- names(which(match.l))
    x.temp <- t(as.matrix(x.temp))
    rownames(x.temp) <- cl.name
  }
  return(x.temp)
}

genMissingAlleleNa <- function(ref.df, allele){
  ref.df$allele <- rep(allele, dim(ref.df)[1])
  ref.df$concordance <- rep(NA, dim(ref.df)[1])
  return(ref.df)
}

matchMeltedAlleles <- function(x.ab.list){
  ab.names <- names(x.ab.list)
  melt.xa <- x.ab.list[[ab.names[1]]]
  melt.xb <- x.ab.list[[ab.names[2]]]

  melt.xa[which(melt.xa$concordance  > 0), 'concordance'] <- 1
  melt.xb[which(melt.xb$concordance  > 0), 'concordance'] <- 1
  all.cl <- c(as.vector(unlist(melt.xa$cell.line)),
              as.vector(unlist(melt.xb$cell.line)))
  all.cl <- unique(all.cl[order(all.cl)])

  for(each.cl in all.cl){
    xa.rows <- which(melt.xa$cell.line %in% each.cl)
    xb.rows <- which(melt.xb$cell.line %in% each.cl)

    if(length(xa.rows) == 0 & length(xb.rows) > 0){
      # If A-allele has no discordancies
      melt.xb[xb.rows, 'match'] <- "A-absent"
    } else if(length(xb.rows) == 0 & length(xa.rows) > 0){
      # If B-allele has no discordancies
      melt.xa[xa.rows, 'match'] <- "B-absent"
    } else {
      xa.vals <- melt.xa[xa.rows, 'concordance']
      xb.vals <- melt.xb[xb.rows, 'concordance']

      xa.vals[is.na(xa.vals)] <- -1
      xb.vals[is.na(xb.vals)] <- -1
      if(!all(xa.vals == xb.vals)){
        # If there is partial match between A- and B-allele discordancies
        melt.xa[xa.rows, 'match'] <- "partial"
        melt.xb[xb.rows, 'match'] <- "partial"
      } else if(all(xa.vals == xb.vals)){
        # If A-allele perfectly matches B-alleles
        melt.xa[xa.rows, 'match'] <- "match"
        melt.xb[xb.rows, 'match'] <- "match"
      } else {
      print(paste("Error", each.cl, sep=""))
      }
    }
  }
  x.ab.list[[ab.names[1]]] <- melt.xa
  x.ab.list[[ab.names[2]]] <- melt.xb
  return(x.ab.list)
}


# Function: addAnnoLines
# Purpose:  Adds annotations and lines to given spoints
# Input:   anno.points <- a named vector of percent identities
# Returns:
addAnnoLines <- function(anno.points, max.height=0.75, x.pos=1, dist.y=0.05, dir='right', segcol=alpha("black", 0.5)){
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
           y1=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y),
           col=segcol)
  # Horizontal line from diagonal to name
  segments(x0=rep(x.pos.1, length(anno.points)),
           y0=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y),
           x1=rep(x.pos.2, length(anno.points)),
           y1=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y),
           col=segcol)
  # Name of Segments
  text(x=rep(x.text, length(anno.points)),
       y=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y),
       labels=anno.names,
       adj=if(dir == 'right') 0 else 1)
}


# Function: matchClAnno
# Purpose:  Given a threshold, it will go through each cell line and determine which cell lines match and which don't match
# Input:  x <- a row in cell.line.anno dataframe, that contains 'cgp.filename' and 'ccle.filename'
#         match_val <- A threshold value (default = 0.8 for snp, 0.6 for cnv)
#         matrix.perc.match <-  a matrix containing all possibel pair-wise concordance scores
# Returns:  List containing matches and mismatches and which cell lines are involved
matchClAnno <- function(x, match_val, matrix.perc.match, fn.headers){
  match.values <- c()

  print(as.character(x['GDSC.cellid']))
  cl.fn <- as.vector(x[fn.headers])
  cl.fn <- cl.fn[!is.na(cl.fn)]
  rm.fn <- ""

  # Checks to make sure all the reported cell.line filenames from fn.headers are in cl.fn
  # removes ones that are not in that list.
  if(length(cl.fn[!(cl.fn %in% colnames(matrix.perc.match))]) > 0){
    rm.fn <- cl.fn[!(cl.fn %in% colnames(matrix.perc.match))]
    #print(paste("Not found: ", cl.fn[!(cl.fn %in% colnames(matrix.perc.match))],
    #            "  --  ", "unique cell id: ", x['unique.cellid'], sep=""))
  }

  # Checks to see if the number of non-matching cl.fn matches the length of cl.fn
  if(length(cl.fn[!(cl.fn %in% colnames(matrix.perc.match))]) == length(cl.fn)){
    #print(paste("None of the dataset filenames were found in the match matrix:\n ", x['unique.cellid'], sep=""))
    match.values <- -1
    cl.fn <- 'NA'

    # Finds all the matching or mismatching based on the annotated matching cl.fn list
  } else {
    if(match_val == "match"){
      cl.fn <- cl.fn[cl.fn %in% colnames(matrix.perc.match)]
      match.values <- matrix.perc.match[cl.fn, cl.fn]
      if((class(match.values) == "numeric") && (length(match.values) == 1)){
        match.values <- as.matrix(match.values)
      }

      # Sets colnames and rownames to dataset names & cell line name for matching annotations
      cl.fn <- as.vector(x[fn.headers])
      # Removes na's and files not in the matrix.perc.match
      fn.hd <- gsub(".filename(.+)?", "", fn.headers[!(cl.fn %in% rm.fn) & !is.na(cl.fn)], perl=TRUE)
      fn.hd <- paste(fn.hd, x['unique.cellid'], sep="_")
      colnames(match.values) <- fn.hd
      rownames(match.values) <- fn.hd

    } else if(match_val == "mismatch"){
      cl.fn <- cl.fn[(cl.fn %in% colnames(matrix.perc.match))]
      cl.mismatch.fn <- colnames(matrix.perc.match)[!(colnames(matrix.perc.match) %in% cl.fn)]
      if(!is.na(table(matrix.perc.match[cl.fn,cl.mismatch.fn] %in% NA)["TRUE"])){
        #print(cl.fn)
      }
      match.values <- matrix.perc.match[cl.fn, cl.mismatch.fn]
    }
  }

  if(class(match.values) == "numeric"){
    match.values <- t(match.values)
    #rownames(match.values) <- cl.fn
  }

  return(as.matrix(match.values))

}

cleanMismatchList <- function(x){
  #Clean up mismatch matrix
  if(!any(dim(x) == 0)){
    x <- x[unique(rownames(x)),unique(colnames(x)), drop=FALSE]

    # Get the index in cell.line.anno for the given ID and the highly correlated Mismatch IDs
    ref.anno.idx <- t(sapply(rownames(x), function(y) which(y == cell.line.anno, arr.ind=TRUE)))
    mismatch.anno.idx <- t(sapply(colnames(x), function(y){
      match.cla <- which(y == cell.line.anno, arr.ind=TRUE)
      if(dim(match.cla)[1] < 1){
        #If the cell.ID is not found in cell.line.anno
        match.cla <- rbind(match.cla, c(NA,NA))
      }
      return(match.cla[1,])
    }))

    # Validate that all ref.annos are the same:
    if(!length(unique(ref.anno.idx[,1])) == 1){
      print(paste("Error: Mismatch IDs do not match for samples: ",
                  paste(unique(rownames(x)), collapse=","), sep=""))
    } else {

      # Flag any mismatches that have different annotations
      ref.cla.row <- unique(ref.anno.idx[,1])
      mm.cla.row <- unique(mismatch.anno.idx[,1])

      if(!all(ref.cla.row == mm.cla.row, na.rm = TRUE)){
        # Subset for just the differently annotated mismatches
        mm.cla.row <- which(!mismatch.anno.idx[,1] == ref.cla.row)
        x <- x[,mm.cla.row, drop=FALSE]

        colnames(x) <- getCellids(mismatch.anno.idx[mm.cla.row,,drop=FALSE])
        rownames(x) <- getCellids(ref.anno.idx)
      } else {
        x <- matrix(ncol=0, nrow=0)
      }
    }
  }
  return(x)
}

getCellids <- function(x){
  fn.headers <- colnames(cell.line.anno)[x[,2]]
  uniq.cellid <- cell.line.anno[x[,1], "unique.cellid"]
  fn.headers <- paste(gsub(".filename(.+)?", "", fn.headers, perl=TRUE), uniq.cellid, sep="_")
  return(fn.headers)
}
}

###############################
#           Main
###############################
#Load and reformat/rename the files
# cell.line.anno$gdsc.filename.x <- gsub(".cel$", ".CEL", cell.line.anno$gdsc.filename.x)
# cell.line.anno$gdsc.filename.y <- gsub(".cel$", ".CEL", cell.line.anno$gdsc.filename.y)
col.index <- which(colnames(cell.line.anno) %in% c('pfizer.filename.y', 'pfizer.filename', 'gdsc.filename.y'))
# colnames(cell.line.anno)[col.index] <- c("pfizer2.filename.y", "pfizer3.filename", "gdsc2.filename.y")

# Match CEL to cell lines and generate match/mismatches
all.headers <- c('cgp.filename', 'ccle.filename', 'gdsc.filename.x', 'pfizer.filename.x',
                 'pfizer2.filename.y', 'pfizer3.filename', 'gdsc2.filename.y')
all.filenames <- c()

# Match CEL to cell lines and generate match/mismatches
low.match.list <- list()
ds.type <- c('cgp.filename', 'ccle.filename', 'gdsc.filename.x', 'pfizer.filename.x',
             'pfizer2.filename.y', 'pfizer3.filename', 'gdsc2.filename.y')

for(ds1 in ds.type){
  low.match.df <- data.frame(matrix(ncol=1, nrow=0))
  low.match.list[[ds1]] <- low.match.df

  for(ds2 in ds.type){
    if(!ds1 == ds2){
      fn.headers <- c(ds1, ds2)  # colnames in cell.line.anno for CGP and CCLE filenames
      pdf(paste(paste(gsub(".filename(.+)?", "", fn.headers, perl=TRUE), sep="", collapse="_"), ".annoplot.pdf", sep=""))
      layout(matrix(c(rep(c(0,3,4,4,2,2,1,0), 9),
                    c(0,0,0,0,5,5,0,0)),
                    ncol=8, nrow=10, byrow=TRUE))
      mat.count <- 1

      for(matrix.perc.match in list(hscr.a2.mat, hscr.a1.mat)){
        colnames(matrix.perc.match) <- gsub("$", ".CEL", colnames(matrix.perc.match))
        rownames(matrix.perc.match) <- gsub("$", ".CEL", rownames(matrix.perc.match))
        print(paste("Cnt: ", ds1, "-", ds2, ": ", mat.count, sep=""))
        na.idx <- which(apply(matrix.perc.match, 1, function(x) all(is.na(x))))
        matrix.perc.match <- matrix.perc.match[-na.idx, -na.idx]

        match.list <- apply(cell.line.anno, 1, function(x) matchClAnno(x, "match", matrix.perc.match, fn.headers))
        names(match.list) <- as.vector(cell.line.anno$unique.cellid)
        mismatch.list <- apply(cell.line.anno, 1, function(x) matchClAnno(x, "mismatch", matrix.perc.match, fn.headers))
        names(mismatch.list) <- as.vector(cell.line.anno$unique.cellid)



        ###############################
        #           Visualization
        ###############################
        #Create match for each cell line density and scatterplot
        match.v <- unlist(match.list, use.names=TRUE)
        match.v <- match.v[which(match.v != 1)]   # cell lines that match to themselves, diagonal
        match.v <- match.v[which(match.v != -1)]  #
        if(length(match.v) > 0){
          match.v <- match.v[-seq(1, length(match.v), by=2)] #Remove duplicated entries
          names(match.v) <- gsub("[1234]$", "", names(match.v), perl=TRUE) #Remove index number of matrix
          dens.mm <- density(unlist(mismatch.list), na.rm=TRUE)
          if(length(match.v) > 2) {
            dens.match <- density(match.v)
          } else {
            dens.match <- density(rnorm(100, mean = 0.9, sd = 0.01))
          }
          #low.match.threshold[data.type] <- quantile(match.v, c(0.05))


          #Create mismatch for each cell line density and scatterplot
          #mismatch.v <- unlist(mismatch.list, use.names=TRUE)
          quant.match.cutoff <- quantile(match.v, 0.05)

          #http://www.r-bloggers.com/absolute-deviation-around-the-median/
          quant.match.cutoff <- median(match.v) - (3*mad(match.v, center=median(match.v), constant=1.4826, na.rm=TRUE))

          if(quant.match.cutoff < low.match.threshold[data.type]) {cutoff.val <- low.match.threshold[data.type] } else { cutoff.val <- quant.match.cutoff}
          #cutoff.val <- 0.7

          if(mat.count == 1){ par(mar=c(0, 0, 4.1, 4.1)) } else { par(mar=c(0, 4.1, 4.1, 0)) }
          #par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
          # Plot Shaded density plots for the match and mismatch distributions
          plot(if(mat.count == 1){ dens.match$y } else { ( dens.match$y * -1 ) },
               dens.match$x,
               type="l", col=match.col,
               ylim=c(if(data.type=='cnv') -1 else 0,1),
               xlab="",
               ylab=if(data.type=='cnv' & mat.count == 2) 'Weighted Pearson Correlation' else if(data.type=='snp' & mat.count == 1) 'Concordance',
               axes=FALSE)
          lines(if(mat.count == 1){ (dens.mm$y * (max(dens.match$y) / max(dens.mm$y))) } else { (dens.mm$y * (max(dens.match$y) / max(dens.mm$y)) * -1) },
                dens.mm$x,
                col=mismatch.col)
          # Add in coloured shades
          with(dens.match, polygon(x=if(mat.count == 1){ dens.match$y } else { ( dens.match$y * -1 ) },
                                   y= dens.match$x,
                                   col=alpha(match.col, 0.5)))
          with(dens.mm, polygon(x=if(mat.count == 1){ (dens.mm$y * (max(dens.match$y) / max(dens.mm$y))) } else { (dens.mm$y * (max(dens.match$y) / max(dens.mm$y)) * -1) },
                                y= dens.mm$x,
                                col=alpha(mismatch.col, 0.5)))

          axis(side = 3, at=c(-100, 100), labels=FALSE) #Top
          axis(side = 1, at=c(-100, 100), labels=FALSE) #Bottom
          if(mat.count == 1){
            axis(side = 4, at=c(-100, 100), labels=FALSE) #solid black bar to the right
          }
          if(mat.count == 2){
            axis(side = 2, at=c(-100, -1, -0.5, 0, 0.5, 1.0, 100), labels=TRUE) # solid black bar + labels to the left
          }

          abline(h = low.match.threshold[data.type], col="grey", lty=5)
          abline(h = quant.match.cutoff, col="dimgrey", lty=3)

          # Plot the individual points for match and mismatch points - highlighting the mismatch ones
          par(mar=c(0, 0, 4.1, 0))
          plot(rep(if(mat.count == 1) 1 else 9, length(unique(match.v))),
               unique(match.v),
               col=alpha(match.col, 0.2),
               ylim=c(if(data.type=='cnv') -1 else 0,1), xlim=c(0,10), las=2,
               axes=FALSE,
               xlab='', ylab="",
               main=paste("hscr-a", if(mat.count == 1) "2" else "1", sep=""))
          #paste(gsub(".filename(.+)?", "", fn.headers, perl=TRUE), sep="", collapse=" - "))
          text(x=if(mat.count == 1) 3.1 else 6.9,
               y=1.0, labels = paste("n = ", length(match.v), sep=""),
               adj=if(mat.count == 1) 0 else 1)

          axis(side = 3, at=c(-100, 100), labels=FALSE) #Top
          axis(side = 1, at=c(-100, 100), labels=FALSE) #Bottom
          if(mat.count == 1){
            axis(side = 2, at=c(-100, 100), labels=FALSE, col="grey") #gray bar to the left
          }
          if(mat.count == 2){
            axis(side = 4, at=c(-100, 100), labels=FALSE, col="grey") # gray bar to the right
          }


          # Add in low-match concordance "matching" points
          low.match <- match.v[which(match.v < cutoff.val)]
          if(length(low.match) < max.num.low.match) max.lowm <- length(low.match) else max.lowm <- max.num.low.match
          low.match <- sort(low.match)[c(0:max.lowm)]
          if(!all(is.na(coi))) low.match <- match.v[which(names(match.v) %in% coi)]

          points(rep(if(mat.count == 1) 1 else 9, length(low.match)),  low.match,
                 pch=if(!all(is.na(coi))) 16 else 1,
                 col=if(!all(is.na(coi))) coi.col  else match.col)
          if(length(low.match) > 0){
            addAnnoLines(low.match,
                         max.height=if(data.type=='cnv') (cutoff.val - 0.05) else 0.75,
                         x.pos=if(mat.count == 1) 1.15 else 8.95,
                         dir=if(mat.count == 1){  'right' } else { 'left' })
          }
          abline(h = low.match.threshold[data.type], col="grey", lty=5)
          abline(h = quant.match.cutoff, col="dimgrey", lty=3)

          if(mat.count == 2){
            par(mar=c(0,0,0,0))
            plot(0, type='n', xlim=c(0,10), ylim=c(0,10), axes=FALSE, ylab='', xlab='')
            legend(0.6,
                   10, # places a legend at the appropriate place
                   c('matching cell lines','mismatching cell lines'), # puts text in the legend
                   lty=c(1,1), # gives the legend appropriate symbols (lines)
                   lwd=c(2.5,2.5),col=c(match.col,mismatch.col)) # gives the legend lines the correct color and width
          }
          all.filenames <- c(all.filenames, as.vector(t(cell.line.anno[which(cell.line.anno$unique.cellid %in% names(low.match)), all.headers])))
          mat.count <- mat.count + 1

          low.match.list[[ds1]] <- t(smartbind(t(low.match.list[[ds1]]), t(data.frame(low.match))))
        } else {
          # If there is no matches found between the two datasets
          low.match <- NA
          low.match.list[[ds1]] <- t(smartbind(t(low.match.list[[ds1]]), t(data.frame(low.match))))
        }
        mismatch.cutoff.list <- lapply(mismatch.list, function(x){
          match.ind <- which(x>as.numeric(quant.match.cutoff), arr.ind=TRUE)
          return(x[match.ind[,1], match.ind[,2], drop=FALSE])
        })
        mismatch.clean.list <-  lapply(mismatch.cutoff.list, cleanMismatchList)

        save(mismatch.clean.list, file=file.path(getwd(), paste("mismatch",
                                                                gsub(".filename.*$", "", ds1),
                                                                gsub(".filename.*$", "", ds2),
                                                                mat.count, "Rdata", sep=".")))
        save(match.list, file=file.path(getwd(), paste("match",
                                                       gsub(".filename.*$", "", ds1),
                                                       gsub(".filename.*$", "", ds2),
                                                       mat.count, "Rdata", sep=".")))
      }
      dev.off()
    } else {
      # if ds2 = ds1
      low.match <- NA
      low.match.list[[ds1]] <- t(smartbind(t(low.match.list[[ds1]]), t(data.frame(low.match)))) #HSCR-A1
      low.match.list[[ds1]] <- t(smartbind(t(low.match.list[[ds1]]), t(data.frame(low.match)))) #HSCR-A2

    }
  }
  low.match.list[[ds1]] <- low.match.list[[ds1]][,-1, drop=FALSE]
  coln.ds <- paste(matrix(rep(gsub(".filename.*$", "", ds.type),2), byrow=TRUE, nrow=2), c("a", "b"), sep="_")
  if(dim(low.match.list[[ds1]])[2] == length(coln.ds)){
    colnames(low.match.list[[ds1]]) <- coln.ds
  } else {
    print("Warning: Could not append colnames based on datasets")
    print(paste(ds.type, collapse=" - "))
    print(dim(low.match.list[[ds1]]))
  }

}
#save(low.match.list, file=file.path(getwd(), "low-match-list.Rdata"))

a1.low.match.list <- lapply(low.match.list, function(x) x[,grep("_a", colnames(x)), drop=FALSE])
a1.low.match.list <- lapply(a1.low.match.list, function(x) rmNaRows(x))
a2.low.match.list <- lapply(low.match.list, function(x) x[,grep("_b", colnames(x)), drop=FALSE])
a2.low.match.list <- lapply(a2.low.match.list, function(x) rmNaRows(x))
names(low.match.list) <- gsub(".filename.*$", "", names(low.match.list))




dir.create(file.path(global.wd, "allele_nc"))
names(low.match.list) <- gsub(".filename.*$", "", names(low.match.list))
lapply(names(low.match.list), function(x) plotNoncorcodantCl(low.match.list[[x]], x,
                                                                file.path(global.wd, "allele_nc")))

# Reports summary statistics based on matching annotation but discordant profiles
if(report.summary){
  all.mm.cl.names <- sapply(low.match.list, function(x) rownames(x))
  # Finds total of unique mismatched cell line pairs in all datasets
  all.mm.cl.v <- levels(factor(as.vector(unlist(all.mm.cl.names))))
  all.mm.cl.v <- all.mm.cl.v[-which(all.mm.cl.v %in% c("T.T", "TT", "KM-H2"))]
  print(paste("Total CL in all datasets with same annotation but discordant profile: ",
              length(all.mm.cl.v), sep=""))

  # Finds total of unique mismatched cell lines pairs per datasets
  cnt.ds.mm <- sapply(all.mm.cl.names, length)
  cnt.ds.filt.mm <- sapply(all.mm.cl.names, function(x) length(x[which(x %in% c("KM-H2", "T.T", "TT"))]))
  print("Dataset level of CL withsame annotation but discordant profile: ")
  print(cnt.ds.mm - cnt.ds.filt.mm)
  pfizer.ds.mm <- all.mm.cl.names[c("pfizer", "pfizer2", "pfizer3")]
  pfizer.cnt.mm <- levels(factor(as.vector(unlist(pfizer.ds.mm))))
  print(paste("Pfizer count: ", length(pfizer.cnt.mm), sep=""))

  print(paste("Total CGP file: ", length(which(gsub(".cel$", "", cell.line.anno$cgp.filename, ignore.case=TRUE)
                                               %in% colnames(hscr.a2.mat))), sep=""))
  print(paste("Total CCLE file: ", length(which(gsub(".cel$", "", cell.line.anno$ccle.filename, ignore.case=TRUE)
                                                %in% colnames(hscr.a2.mat))), sep=""))
  print(paste("Total GDSC file: ", length(which(gsub(".cel$", "", cell.line.anno$gdsc.filename, ignore.case=TRUE)
                                                %in% colnames(hscr.a2.mat))), sep=""))
  print(paste("Total Pfizer file: ", length(which(gsub(".cel$", "", cell.line.anno$pfizer.filename, ignore.case=TRUE)
                                                  %in% colnames(hscr.a2.mat))), sep=""))
}

# Reports a summary of just cell line name and dataset to be used in phenotype plotting (script unnamed as of 02/11/2016)
ds.corr.simple <- lapply(low.match.list, function(ds.corr){
  ds.disc <- apply(ds.corr, 1, function(row.corr){
    non.na.idx <- which(!is.na(row.corr))
    ds.disc <- gsub("_[ab]", "", colnames(ds.corr)[non.na.idx])
    ds.disc <- unique(gsub("[123]$", "", ds.disc))
    return(ds.disc)
  })
  ds.id <- list()
  for(each.cl in c(1:length(ds.disc))){
    ds.id[[rownames(ds.corr)[each.cl]]] <- sapply(ds.disc[each.cl], function(x) {
                                                    return(paste(x, rownames(ds.corr)[each.cl], sep="_"))
                                                  })
  }
  return(ds.id)
})
#save(ds.corr.simple, file=file.path("rdata", "low_match_list.simple.Rdata"))
