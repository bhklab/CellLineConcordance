library(AnnotationDbi)
library(org.Hs.eg.db)
library(plyr)
library(vioplot)
library(PharmacoGx)
library(scales)
library(gtools)
library(intervals)
library(preprocessCore)
library(Hmisc)
library(ggplot2)

## Load Cell line annotation matrix
global.curation.dir <- '~/git/cnv_fingerprint/v2.0/data';
out.rdata <- 'merged_annotations.Rdata';
load(file.path(global.curation.dir, out.rdata))

## Load CNV correlation matrix
cnv.matrix.dir <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/cnv_analysis'
cnv.matrix.file <- file.path(cnv.matrix.dir, "input", "hscra1.weightedCorrDf.Rdata")
load(cnv.matrix.file)
hscr.a1.mat <- weighted.corr.df
cnv.matrix.file <- file.path(cnv.matrix.dir, "input", "hscra2.weightedCorrDf.Rdata")
load(cnv.matrix.file)
hscr.a2.mat <- weighted.corr.df

## Load Drug sensitivity data
GDSC <- downloadPSet("GDSC")
CCLE <- downloadPSet("CCLE")
commonClAuc <- intersectPSet(list('CCLE'=CCLE,
                                  'GDSC'=GDSC), 
                             intersectOn=c("cell.lines", "drugs"))
drug.auc.list <- list()
drug.auc.list[['GDSC']] <- summarizeSensitivityProfiles(
  pSet=commonClAuc$GDSC,
  sensitivity.measure='auc_recomputed',
  summary.stat='median')
drug.auc.list[['CCLE']] <- summarizeSensitivityProfiles(
  pSet=commonClAuc$CCLE,
  sensitivity.measure='auc_recomputed',
  summary.stat='median')


################
#### FUNCTIONS
## Multiplot, plot multiple ggplot2: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, rows=NA, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  if(is.list(...)) {
    plots <- c(..., plotlist)
  } else {
    plots <- c(list(...), plotlist)
  }
  
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, 
                     nrow = if(is.na(rows)) ceiling(numPlots/cols) else rows)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#### Compute the ECDFs for Deltas in drug AUC between cell lines
getDeltaAuc <- function(auc.x, auc.y, cell.ids, analysis.type='concordant', filt=FALSE){
  if(nrow(auc.x) == nrow(auc.y)){
    if(analysis.type == 'concordant'){
      delta.auc.list <- list()
      auc.x <- auc.x[,which(colnames(auc.x) %in% cell.ids), drop=FALSE]
      auc.y <- auc.y[,which(colnames(auc.y) %in% cell.ids), drop=FALSE]
      
      for(each.drug in rownames(auc.x)){
        combine.drug <- rbind(auc.x[each.drug, ,drop=FALSE], 
                              auc.y[each.drug, ,drop=FALSE])
        if(filt) combine.drug[combine.drug < 0.20] <- NA
        delta.auc.tmp <- abs(apply(combine.drug, 2, diff))
        delta.auc.list[[each.drug]] <- delta.auc.tmp
      }
      delta.auc <- do.call("rbind", delta.auc.list)
    } else if(analysis.type=="discordant"){
      # Generate correlation of all expressions between Dataset X and Dataset Y
      # for cell lines that are annotated differently. 
      print("Obtaining difference of AUC for differently annotated cell lines")
      if(is.na(cell.ids)){
        delta.auc <- sapply(colnames(auc.x), function(each.cl.x){
          combine.expr.row <- sapply(colnames(auc.y), function(each.cl.y){
            if(each.cl.x != each.cl.y){
              combine.drug <- rbind(auc.x[,each.cl.x,drop=FALSE], 
                                    auc.y[,each.cl.y,drop=FALSE])
              delta.auc.tmp <- abs(apply(combine.drug, 2, diff))
            } else {
              delta.auc.tmp <- NA
            }
          })
          return(combine.expr.row)
        })
      } else {
        print("Error: Intersect based on common cell lines and drugs first")
        delta.auc <- NA
      }
    }
  }
  return(delta.auc)
}

# Function: matchClAnno
# Purpose:  Given a threshold, it will go through each cell line and determine which cell lines match and which don't match
# Input:  x <- a row in cell.line.anno dataframe, that contains 'cgp.filename' and 'ccle.filename'
#         match_val <- A threshold value (default = 0.8 for snp, 0.6 for cnv)
#         matrix.perc.match <-  a matrix containing all possibel pair-wise concordance scores
# Returns:  List containing matches and mismatches and which cell lines are involved
matchClAnno <- function(x, match_val, matrix.perc.match, fn.headers){
  match.values <- c()
  
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

getCorrVals <- function(hscr_list, fn.headers, cell.line.anno){
  corr_vals <- list()
  m_cnt <- 1
  for(matrix.perc.match in hscr_list){
    colnames(matrix.perc.match) <- gsub("$", ".CEL", colnames(matrix.perc.match))
    rownames(matrix.perc.match) <- gsub("$", ".CEL", rownames(matrix.perc.match))
    match.list <- apply(cell.line.anno, 1, function(x) matchClAnno(x, "match", matrix.perc.match, fn.headers))
    names(match.list) <- as.vector(cell.line.anno$unique.cellid)
  
    ##
    #Create match for each cell line density and scatterplot
    match.v <- unlist(match.list, use.names=TRUE)
    match.v <- match.v[which(match.v != 1)]   # cell lines that match to themselves, diagonal
    match.v <- match.v[which(match.v != -1)]  #
    if(length(match.v) > 0){
      match.v <- match.v[-seq(1, length(match.v), by=2)] #Remove duplicated entries
      names(match.v) <- gsub("[1234]$", "", names(match.v), perl=TRUE) #Remove index number of matrix
    }
    corr_vals[[m_cnt]] <- match.v
    m_cnt <- m_cnt + 1
  }
  return(corr_vals)
}


################
#### MAIN
ds1 <- 'gdsc.filename.x'
ds2 <- 'ccle.filename'
fn.headers <- c(ds1, ds2)  # colnames in cell.line.anno for CGP and CCLE filenames
corr_vals <- getCorrVals(hscr_list=list(hscr.a2.mat, hscr.a1.mat), fn.headers, cell.line.anno)
corr_vals_df <- do.call("rbind", corr_vals)

# Delta AUC for matching cell lines that have discordant CNVs
cnvDiff.delta.auc <- getDeltaAuc(auc.x=drug.auc.list[['CCLE']], 
                                 auc.y=drug.auc.list[['GDSC']], 
                                 cell.ids=names(corr_vals[[1]]), analysis.type='concordant',
                                 filt=FALSE)
aucCnv_df <- rbind.fill(as.data.frame(cnvDiff.delta.auc), as.data.frame(corr_vals_df))
rownames(aucCnv_df) <- c(rownames(cnvDiff.delta.auc), "cnv1", "cnv2")


plist <- list()
for(each_drug in rownames(cnvDiff.delta.auc)){
  each_drug_filt <- t(aucCnv_df[c(each_drug, "cnv1", "cnv2"), ,drop=FALSE])
  each_drug_filt <- as.data.frame(each_drug_filt[which(!is.na(each_drug_filt[,1])),])
  colnames(each_drug_filt) <- c("deltaAuc", "cnv1", "cnv2")
  each_drug_filt$avgcnv <- unlist(apply(each_drug_filt[,c("cnv1", "cnv2")], 1, mean))
  each_drug_filt$Bins <- 0
  each_drug_filt[with(each_drug_filt, 
                      which(avgcnv < 0.8 | deltaAuc > quantile(deltaAuc, 0.9))),'Bins'] <- 1
  # plot(deltaAuc ~ avgcnv, type='n', data=each_drug_filt, col="black", main=each_drug,
  #      ylim=c(0,1), xlim=c(0.3, 1))
  # text(deltaAuc ~ avgcnv, data=each_drug_filt, labels=rownames(each_drug_filt), col="black")
  # text(x=0.3, y=0.1, adj=0, labels=paste("r =", round(cor(each_drug_filt)[1,4], 3)))
  
  plist[[each_drug]] <- ggplot(each_drug_filt, aes(x=avgcnv, y=deltaAuc, colour=Bins)) +
    geom_point(alpha = 0.5 ) +
    xlim(0.3, 1) +
    ylim(0, 1) +
    labs(title = each_drug, x = "", y = "") +
    theme(legend.position="none")  
  
}
pdf("drugAuc-Cnv.pdf")
multiplot(plist, cols=4)
dev.off()
cat(file.path(getwd(), "drugAuc-Cnv.pdf\n"))