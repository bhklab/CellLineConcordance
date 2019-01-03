library(gplots)
library(gtools)
library(scales)
library(reshape)
library(ggplot2)
library("RColorBrewer")

plotNoncordantCl <- function(x, title, target, binarize=FALSE, 
                             nonmatch.threshold=0.7, ambiguous.threshold=0.8,
                             order=TRUE){
  myPalette <- colorRampPalette(c(mismatch.col, "white"), bias=1, interpolate="linear")
  
  x.m <- melt(as.matrix(x))
  colnames(x.m) <- c("cell.line", "dataset", "concordance")
  x.m$dataset <- factor(x.m$dataset, levels=target)
  if(order) x.m$cell.line <- factor(x.m$cell.line, levels = rownames(x))
  
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

orderRows <- function(dat, order.metric='min', ret.na=FALSE, row.idx=NULL){
  print(paste0("Ordering metric: ", order.metric))
  na.idx <- which(apply(dat, 1, function(x) all(is.na(x))))
  if(length(na.idx) > 0){
    naomit.dat <- dat[-na.idx,]
    onlyna.dat <- dat[na.idx,]
  } else {
    naomit.dat <- dat
    onlyna.dat <- NULL
  }
  
  
  
  switch(order.metric,
         min={
           ord.idx <- order(apply(naomit.dat, 1, min, na.rm=TRUE), decreasing = TRUE)
         },
         mean={
           ord.idx <- order(apply(naomit.dat, 1, mean, na.rm=TRUE), decreasing = TRUE)
         },
         max={
           ord.idx <- order(apply(naomit.dat, 1, max, na.rm=TRUE), decreasing = TRUE)
         },
         custom={
           ord.idx <- match(row.idx, rownames(naomit.dat))
         }
         )
  
  dat <- naomit.dat[ord.idx,]
  if(ret.na & (length(na.idx) > 0)) dat <- rbind(dat, onlyna.dat)
  dat
}

genSummMetrics <- function(){
  load(paste0('nBraw', ".GNE_matchdf.Rdata"))
  nb.anno.df <- orderRows(match.anno.df, order.metric='min', ret.na=TRUE)
  load(paste0('nAraw', ".GNE_matchdf.Rdata"))
  na.anno.df <-  orderRows(match.anno.df, order.metric='min', ret.na=TRUE)
  
  nAB.list <- list(as.numeric(na.anno.df), as.numeric(nb.anno.df))
  means <- round(sapply(nAB.list, mean, na.rm=TRUE),2)
  sds <- round(sapply(nAB.list, sd, na.rm=TRUE),2)
  print(paste("Pearson correlation of: ",  means,  "+/-", sds))
  
  threshold <- function(x, case, thresh=0.8, thresh.lo=0.6) {
    if(all(is.na(x))){
      res <- FALSE
    } else {
      switch(case,
             "ge"= res <- all(x >= thresh, na.rm=TRUE),
             "le" = res <- all(x <= thresh, na.rm=TRUE),
             "gal" = res <- all(x < thresh & x > thresh.lo, na.rm=TRUE))
    }
    return(res)
  }
  match.a <- table(apply(na.anno.df, 1, function(x) threshold(x, 'ge', 0.8)))[2]
  match.b <- table(apply(nb.anno.df, 1, function(x) threshold(x, 'ge', 0.8)))[2]
  
  ambig.a <- table(apply(na.anno.df, 1, function(x) threshold(x, 'gal', 0.8, 0.6)))[2]
  ambig.b <- table(apply(nb.anno.df, 1, function(x) threshold(x, 'gal', 0.8, 0.6)))[2]
  
  disc.a <- table(apply(na.anno.df, 1, function(x) threshold(x, 'le', 0.6)))[2]
  disc.b <- table(apply(nb.anno.df, 1, function(x) threshold(x, 'le', 0.6)))[2]
  
  na.a <- table(apply(na.anno.df, 1, function(x) all(is.na(x))))[2]
  na.b <- table(apply(nb.anno.df, 1, function(x) all(is.na(x))))[2]
  
  summ.df<- data.frame("A" = c(match.a, ambig.a, disc.a, na.a),
                       "B" = c(match.b, ambig.b, disc.b, na.b))
  summ.df
}

mismatch.col <- 'magenta'
match.col <- 'brown'


#### COPY NUMBER
setwd("/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/rdataFiles/cnv/output")
out.dir <- "/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/rdataFiles/cnv/output/GNE_summary_plots"

load(paste0('nBraw', ".GNE_matchdf.Rdata"))
match.anno.df <-  orderRows(match.anno.df, order.metric='min', ret.na=TRUE)
names.x <- rownames(match.anno.df)


gen.summ.metrics <- FALSE
if(gen.summ.metrics) genSummMetrics()
        

for (hscr in c("nAraw", "nBraw")) {
  #hscr <- 'nAraw'
  #load(paste0(hscr, ".GNE.corr.Rdata"))
  load(paste0(hscr, ".GNE_matchdf.Rdata"))
  
  #match.anno.df <-  orderRows(match.anno.df, order.metric='min', ret.na=TRUE)
  match.anno.df <-  orderRows(match.anno.df, order.metric='custom', ret.na=FALSE, row.idx=names.x)
  pdf(file.path(out.dir, paste0(hscr, ".conc.pdf")), height = 100)
  colnames(match.anno.df) <- toupper(colnames(match.anno.df))
  plotNoncordantCl(match.anno.df, paste0("GNE-", hscr), colnames(match.anno.df),
                   binarize=TRUE, 0.65, 0.8, TRUE)
  plotNoncordantCl(match.anno.df, paste0("GNE-", hscr), colnames(match.anno.df),
                   binarize=FALSE, 0.65, 0.8, TRUE)
  dev.off()
  
  cat(paste0("scp quever@mordor:", file.path(out.dir, paste0(hscr, ".conc.pdf .\n"))))
  
}



#### GENOTYPE
setwd("/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2017_GNE/snp_fingerprint/output")
out.dir <- getwd()
load("GNE.genotypeConcordance.Rdata")

#match.anno.df <-  orderRows(match.anno.df, order.metric='min', ret.na=TRUE)
match.anno.df <-  orderRows(match.anno.df, order.metric='custom', ret.na=FALSE, row.idx=names.x)

pdf(file.path(out.dir, "genotype.conc.pdf"), height = 100)
colnames(match.anno.df) <- toupper(colnames(match.anno.df))
plotNoncordantCl(match.anno.df, "GNE-genotype", colnames(match.anno.df),
                 binarize=TRUE, 0.7, 0.8, TRUE)
plotNoncordantCl(match.anno.df, "GNE-genotype", colnames(match.anno.df),
                 binarize=FALSE, 0.7, 0.8, order=TRUE)
dev.off()
cat(paste0("scp quever@mordor:", file.path(out.dir, "genotype.conc.pdf"), " .\n"))

