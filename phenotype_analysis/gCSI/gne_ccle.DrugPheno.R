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
library(vioplot)

load("/Users/rquevedo/git/reference/l1000/l1000.Rdata")

####################
## FUNCTIONS
getTissueCommonCl <- function(ds1, ds2, targ.cl){
  getTissue <- function(targ.cl, ds){
    row.idx <- grep(targ.cl, rownames(ds@curation$tissue))
    ds@curation$tissue[row.idx,]$unique.tissueid
  }
  getCellIds <- function(tissue.id, ds){
    tissue.idx <- grep(tissue.id, 
                       ds@curation$tissue$unique.tissueid)
    rownames(ds@curation$tissue)[tissue.idx]
  }
  cell.tissue.1 <- getCellIds(getTissue(targ.cl, ds1), ds1)
  cell.tissue.2 <- getCellIds(getTissue(targ.cl, ds2), ds2)
  intersect(cell.tissue.1, cell.tissue.2)
}

# https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

summKSstat <- function(summ.aac, Y.drug.aac){
  kstat <- pval <- c();
  
  for (drug.id in names(Y.drug.aac)) {
    seq1 <- summ.aac[[drug.id]][[1]]
    seq2 <- summ.aac[[drug.id]][[2]]
    ks <- ks.test(seq1, seq2)
    pval <- c(pval, ks$p.value)
    kstat <- c(kstat, ks$statistic)
  }
  
  
  return(list("D"=kstat,
              "p"=pval))
}

getXVals <- function(refmat, ids.x){
  values.x <- sapply(ids.x, function(x) grep(paste0("^", x, "$"), names(refmat)))
  values.x <- refmat[values.x]
  values.x
}

#' summMatchMatrix
#'
#' @param matx A matrix of concordances/correlations/associations between two datasets.  The diagonal should be the matching cell lines between datasets 
#' @param mat.to.vector Boolean to turn all matrix values into a vector (default=TRUE)
#' @param return.style Returning the "Match" and "Nonmatching" as either a 'data.frame', 'list', or 'melt'
#'
#' @return
#' @export
#'
#' @examples
summMatchMatrix <- function(matx, mat.to.vector=TRUE, return.style='list'){
  matching.diag <- diag(matx)
  diag(matx) <- NA
  if(mat.to.vector) matx <- as.numeric(matx)
  summ.metric <- list("Matching"=matching.diag,
                      "Nonmatching"=matx)
  if(mat.to.vector){
    if(return.style == 'data.frame' | return.style == 'melt'){
      dat <- data.frame("Matching"=c(summ.metric[[1]], 
                                     rep(NA, length(summ.metric[[2]]) - 
                                           length(summ.metric[[1]]))),
                        "Nonmatching"=summ.metric[[2]])
      if(return.style == 'melt') {
        dat <- melt(dat)
        dat <- dat[-which(is.na(dat$value)),]
      }
      # ggplot(dat,aes(x=variable,y=value)) +
      #   geom_violin() +
      #   geom_text(aes(y=max(value,na.rm=TRUE)/2,label='test'))
    } else {
      dat <- summ.metric
    }
  }
  
  return(dat)
}

subsetExprMat <- function(expr.mat, expr.ids, input.type='ENTREZID', 
                          expr.idx=NULL, rm.na=TRUE){
  if(is.null(expr.idx)){
    print("Using org.Hs.eg.db to subset expression matrix...")
    conv.ens.id <- mapIds(org.Hs.eg.db,
                          keys=as.character(expr.ids),
                          column="ENSEMBL",
                          keytype=input.type,
                          multiVals="first")
    expr.idx <- match(paste(conv.ens.id, "at", sep="_"), rownames(expr.mat))
    expr.mat <- expr.mat[expr.idx,]
  } else {
    expr.mat <- expr.mat[intersect(expr.idx, rownames(expr.mat)),]
  }
  
  
  if(rm.na){
    na.idx <- which(apply(expr.mat, 1, function(x) all(is.na(x))))
    if(length(na.idx) > 0) expr.mat <- expr.mat[-na.idx,]
  }
  return(expr.mat)
}
## PharmacoGX Default Wrappers
getExprVals <- function(intersect.list, dataset=NULL, data.type='rna', genes=NA){
  if(!is.null(dataset)) ds.vals <- intersect.list[[dataset]] else ds.vals <- intersect.list
  if(all(is.na(genes))){
    summarizeMolecularProfiles(ds.vals,
                               cellNames(ds.vals),
                               mDataType=data.type,
                               verbose=FALSE)
  } else {
    summarizeMolecularProfiles(ds.vals,
                               cellNames(ds.vals),
                               mDataType=data.type,
                               genes,
                               verbose=FALSE)
  }
}

getDrugVals <- function(intersect.list, dataset=NULL, data.type='auc_recomputed'){
  if(!is.null(dataset)) ds.vals <- intersect.list[[dataset]] else ds.vals <- intersect.list
  summarizeSensitivityProfiles(pSet=ds.vals,
                               sensitivity.measure=data.type,
                               summary.stat='median')
}

plotBayesFactor <- function(seq1, seq2, values.x, 
                            f1col='blue', f2col='red', 
                            kde.plot=TRUE, 
                            annotate.plot=TRUE,
                            dot.plot=FALSE,
                            upp.y=2, low.y=-0.2,
                            lty=2, alpha.val=0.6,
                            cex.range=c(2, 5),
                            ramp.cols=list("1" = "black", "-1" = "red"),
                            ...){
  ks.pval <- ks.test(seq1, seq2)$p.val
  
  # Generates the KDE for the match/nonmatch curves
  if(kde.plot){
    dseq1 <- density(seq1, na.rm=TRUE)
    dseq2 <- density(seq2, na.rm=TRUE)
    
    f1 <- approxfun(dseq1)
    f2 <- approxfun(dseq2)
    
    max.p <- max(c(dseq1$y, dseq2$y), na.rm=TRUE)
    scale.val <- upp.y / max.p
    
    plot(0, type='n', xlim=c(-0.2,1), ylim=c(low.y, upp.y),
         ylab='Density', yaxt='n', ...)
    axis(side = 2, at = seq(0, upp.y, by=upp.y/5), labels=round(seq(0, max.p, by=max.p/5), 2), las=2)
    lines(dseq1$x, dseq1$y * scale.val, col=alpha(f1col, alpha.val), lwd=2, yaxt='n', ...)
    lines(dseq2$x, dseq2$y * scale.val, col=alpha(f2col, alpha.val), lwd=2, yaxt='n', ...)
    text(x=rep(-0.2, 2), y=c(upp.y, low.y), adj=0,
         labels=c("Match", "Mismatch"), col=c(f1col, f2col))
    
    #KS-Statistic
    text(x=0.77, y=(upp.y-0.1), labels="KS-Statistic:", adj=1, cex=0.9)
    text(x=0.8, y=(upp.y-0.1), labels=paste0("p= ", round(ks.pval,4)), adj=0, cex=0.9)
  }

  # Calculate the likelihood ratio and add it to existing plot
  bf.vals <- list()
  if(annotate.plot){
    row.cnt <- 1
    for(each.col in seq_along(values.x)){
      # Obtain post estimates for probability of x in both curves
      x <- values.x[each.col]
      p2 <- f2(x)
      p1 <- f1(x)
      if(!is.na(x) && is.na(p2)) p2 <- 0
      if(!is.na(x) && is.na(p1)) p1 <- 0
      bf <- log(p2/p1)
      bf.id <- names(values.x)[each.col]
      bf.vals[[bf.id]] <- bf
      
      # Plot the point prob estimates and Likelihood Ratios
      if(!kde.plot || is.na(x)) next
      text(x=0.77, y=(upp.y - 0.1 - row.cnt*0.1), adj=1,
           labels=paste0(bf.id, ":"), cex=0.9)
      text(x=0.8, y=(upp.y - 0.1 - row.cnt*0.1), adj=0, cex=0.9,
           labels=paste0("LR= ", if(is.na(x)) 'NA' else round(bf, 2)))
      
      if(is.na(x)) next
      
      y.pos <-  c((p1 * scale.val), (p2 * scale.val))
      points(rep(x, 2), y.pos, col=c(f1col, f2col), pch=19)
      rect(xleft = c(x-0.01, x), ybottom = c(0,0),
           xright = c(x, x+0.01), ytop = y.pos,
           col=alpha(c(f1col, f2col), 0.5), border = FALSE)
      text(x=rep(x, 2), y=c(upp.y, low.y), col=c(f1col, f2col),
           labels=c(round(p1, 2), round(p2, 2)), cex=0.9)
      # text(x=rep(x, 2), y=c(upp.y, low.y), col=c(f1col, f2col),
      #      labels=names(x), cex=0.9)
      row.cnt <- row.cnt + 1
    }
  }
  
  # Bayes factor for each set of points
  if(dot.plot){
    plot(0, type='n', axes=FALSE, ylim=c(0, length(values.x) + 0.5), xlim=c(-0.2,1), ylab="", xlab="")
    axis(side = 1)
    axis(side=2, at=c(1:length(values.x)), labels = names(values.x), las=2)
    for(each.col in seq_along(values.x)){
      # Obtain post estimates for probability of x in both curves
      x <- values.x[each.col]
      dseq1 <- density(seq1, na.rm=TRUE)
      dseq2 <- density(seq2, na.rm=TRUE)
      
      f1 <- approxfun(dseq1)
      f2 <- approxfun(dseq2)
      
      
      p2 <- f2(x)
      p1 <- f1(x)
      if(!is.na(x) && is.na(p2)) p2 <- 0
      if(!is.na(x) && is.na(p1)) p1 <- 0
      bf <- round(log(p2/p1), 1)
      
      if(!is.na(bf)){
        if(!is.list(ramp.cols)){
          ramp.cols <- list("1" = ramp.cols[1],
                            "-1" = ramp.cols[2])
        }
        cex.val <- abs(bf)
        if(cex.val >= max(cex.range)) cex.val <- max(cex.range)
        if(cex.val <= min(cex.range)) cex.val <- min(cex.range)
        
        points(x, each.col, pch=16, cex=cex.val, col=ramp.cols[[as.character(sign(bf))]])
      }
    }
  }

  return(list('ks'=ks.pval,
              'bf'=do.call("cbind", bf.vals)))
}


####################
## MAIN
#### Obtain all PharmacoGX Data
#Download PSets
availablePSets()
gCSI <- downloadPSet("gCSI")
CTRPv2 <- downloadPSet("CTRPv2")
GDSC <- downloadPSet("GDSC")
GDSC1000 <- downloadPSet("GDSC1000")
#load("~/Desktop/bhk_lab/data/pharmacogx/GDSC.RData")
CCLE <-downloadPSet("CCLE")




CTRPv2@curation$tissue <- CTRPv2@curation$cell  # Temporary fix
commonClAuc <- intersectPSet(list('CCLE'=CCLE,
                                  'gCSI'=gCSI), 
                             intersectOn=c("cell.lines", "drugs"))
commonGenes <- intersect(fNames(gCSI, "rnaseq"),
                         fNames(CCLE,"rna"))
commonCl <- intersectPSet(list('CCLE'=CCLE,
                               'gCSI'=gCSI),  
                          intersectOn=c("cell.lines"))

tissueCls <- intersect(cellNames(commonClAuc[[1]]),
                       getTissueCommonCl(gCSI, CCLE, 'PSN1'))
commonTissueCl <- intersectPSet(list('CCLE'=CCLE,
                                     'gCSI'=gCSI),
                                intersectOn = c("drugs", "cell.lines"),
                                cells = tissueIds)

# Find the Expression values for overlapping cell lines:
expr.list <- list()
commonTissueCl
expr.list[['gCSI.tissue']] <- getExprVals(commonTissueCl, 'gCSI', "rnaseq", commonGenes)
expr.list[['CCLE.tissue']] <- getExprVals(commonTissueCl, 'CCLE', "rna", commonGenes)
expr.list[['gCSI']] <- getExprVals(commonCl, 'gCSI', "rnaseq", commonGenes)
expr.list[['CCLE']] <- getExprVals(commonCl, 'CCLE', "rna", commonGenes)
expr.list[['CCLE.all']] <- getExprVals(CCLE, NULL, "rna", commonGenes)
expr.list[['gCSI.all']] <- getExprVals(gCSI, NULL, "rnaseq", commonGenes)

# Find the L1000 Genes Expression values for overlapping cell lines:
expr.sub.list <- list()
expr.sub.list[['gCSI']] <- subsetExprMat(exprs(expr.list[['gCSI']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)
expr.sub.list[['CCLE']] <- subsetExprMat(exprs(expr.list[['CCLE']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)
expr.sub.list[['gCSI.tissue']] <- subsetExprMat(exprs(expr.list[['gCSI.tissue']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)
expr.sub.list[['CCLE.tissue']] <- subsetExprMat(exprs(expr.list[['CCLE.tissue']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)

# Find the Drug Recomputed AUC values for overlapping cell lines:
drug.auc.list <- list()
drug.auc.list[['gCSI']] <- getDrugVals(commonClAuc, 'gCSI')
drug.auc.list[['CCLE']] <- getDrugVals(commonClAuc, 'CCLE')
drug.auc.list[['gCSI.tissue']] <- getDrugVals(commonTissueCl, 'gCSI')
drug.auc.list[['CCLE.tissue']] <- getDrugVals(commonTissueCl, 'CCLE')
drug.auc.list[['CCLE.all']] <- getDrugVals(CCLE)
drug.auc.list[['gCSI.all']] <- getDrugVals(gCSI)


### Generating the corr/deltaAAC matrices
# L1000 Gene expression
Y.l1000 <- cor(as.matrix(expr.sub.list[['gCSI']]), 
               as.matrix(expr.sub.list[['CCLE']]), 
               use = "pairwise.complete.obs", method = "pearson")
Y.l1000 <- cor(as.matrix(expr.sub.list[['gCSI.tissue']]), 
               as.matrix(expr.sub.list[['CCLE.tissue']]), 
               use = "pairwise.complete.obs", method = "pearson")

summ.expr.l1000 <- summMatchMatrix(Y.l1000, TRUE, 'list')

# All gene expression
Y <- cor(as.matrix(exprs(expr.list[['gCSI.tissue']])), 
         as.matrix(exprs(expr.list[['CCLE.tissue']])), 
         use = "pairwise.complete.obs", method = "pearson")
summ.expr <- summMatchMatrix(Y, TRUE, 'list')

# All drug recomputed delta-AAC
Y.drug.aac <- lapply(seq_along(rownames(drug.auc.list[['gCSI.tissue']])), function(each.row){
  Y.tmp <- apply(drug.auc.list[['CCLE.tissue']][each.row, ,drop=FALSE], 2, function(each.col){
    abs(each.col - drug.auc.list[['gCSI.tissue']][each.row, ])
  })
  rownames(Y.tmp) <- colnames(drug.auc.list[['gCSI.tissue']])
  colnames(Y.tmp) <- colnames(drug.auc.list[['CCLE.tissue']])
  Y.tmp
})
names(Y.drug.aac) <- rownames(drug.auc.list[['gCSI.tissue']])
summ.aac <- lapply(Y.drug.aac, summMatchMatrix, mat.to.vector=TRUE, return.style='list')



### Plotting of Likelihood Ratio
out.id <- 'tissue'

switch(out.id,
       tissue={
         ids.x <- 'PSN1'
       },
       cnvDisc={
         ids.x <- c("SR", "PSN1", "NCI-H2052")
         #ids.x <- c("PSN1")
       },
       cnvAmbig={
         ids.x <- c("MOLP-8", "NALM-6", "KYSE-510")
       },
       cnvMatch={
         ids.x <- c("NUGC-3", "SW620", "KARPAS-620")
       })

pdf(paste0("~/Desktop/bhk_lab/results/GNE/", out.id, ".pdf"))
values.x <- getXVals(summ.expr.l1000[[1]], ids.x)
boxplot(summ.expr.l1000)
plotBayesFactor(summ.expr.l1000[[1]], summ.expr.l1000[[2]],
                values.x = values.x,
                xlab='Delta AAC', main="L1000 Gene expression")

values.x <- getXVals(summ.expr[[1]], ids.x)
boxplot(summ.expr)
plotBayesFactor(summ.expr[[1]], summ.expr[[2]],
                values.x = values.x,
                xlab='Correlation', main="Gene expression")

ks.summ <- summKSstat(summ.aac, Y.drug.aac)
ks.summ <- sapply(ks.summ, function(x) c("mean"=mean(x), "sd"=sd(x)))

for (drug.id in names(Y.drug.aac)) {
  values.x <- getXVals(summ.aac[[drug.id]][[1]], ids.x)

  
  boxplot(summ.aac[[drug.id]])
  plotBayesFactor(summ.aac[[drug.id]][[1]], summ.aac[[drug.id]][[2]],
                  values.x = values.x,
                  xlab='Delta AAC', main=drug.id)
}
dev.off()



### Plotting for paper figure
out.id <- 'all_LR_fig'
pdf(paste0("~/Desktop/bhk_lab/results/GNE/", out.id, "2.pdf"), width = 14, height = 14)
split.screen(c(4,7))

# Plot KDE
setScreen <- function(screen.idx, mar=c(2,2,1,1)){
  screen(screen.idx)
  par(mar=mar)
  screen.idx <- screen.idx + 1
}
out.idx <- 0
screen.st.col <- (7 * (out.idx)) + 1
screen.st.col <- setScreen(screen.st.col)
plotBayesFactor(summ.expr.l1000[[1]], summ.expr.l1000[[2]],
                xlab='Delta AAC', main="L1000 Gene expression",
                lty=1, alpha.val=1, f1col="black",
                kde.plot=TRUE, annotate.plot = FALSE, dot.plot = FALSE)

for (drug.id in names(Y.drug.aac)) {
  screen.st.col <- setScreen(screen.st.col)
  plotBayesFactor(summ.aac[[drug.id]][[1]], summ.aac[[drug.id]][[2]],
                  xlab='Delta AAC', main=drug.id,
                  lty=1, alpha.val=1, f1col="black",
                  kde.plot=TRUE, annotate.plot = FALSE, dot.plot = FALSE)
}


# Plot each individual dot plot
all.out.ids <- c("cnvDisc", "cnvAmbig", "cnvMatch")
for (out.idx in seq_along(all.out.ids)) {
  screen.st.col <- (7 * (out.idx)) + 1
  
  out.id <- all.out.ids[out.idx]
  switch(out.id,
         cnvDisc={
           ids.x <- c("SR", "PSN1", "NCI-H2052")
           #ids.x <- c("PSN1")
         },
         cnvAmbig={
           ids.x <- c("MOLP-8", "NALM-6", "KYSE-510")
         },
         cnvMatch={
           ids.x <- c("NUGC-3", "SW620", "KARPAS-620")
         })
  
  values.x <- getXVals(summ.expr.l1000[[1]], ids.x)
  screen.st.col <- setScreen(screen.st.col)
  plotBayesFactor(summ.expr.l1000[[1]], summ.expr.l1000[[2]],
                  values.x = rev(values.x),
                  xlab='Delta AAC', main="L1000 Gene expression",
                  lty=1, alpha.val=1, f1col="black",
                  kde.plot=FALSE, annotate.plot = FALSE, dot.plot = TRUE)
  
  for (drug.id in names(Y.drug.aac)) {
    values.x <- getXVals(summ.aac[[drug.id]][[1]], ids.x)
    
    screen.st.col <- setScreen(screen.st.col)
    plotBayesFactor(summ.aac[[drug.id]][[1]], summ.aac[[drug.id]][[2]],
                    values.x = rev(values.x),
                    xlab='Delta AAC', main=drug.id,
                    lty=1, alpha.val=1, f1col="black",
                    kde.plot=FALSE, annotate.plot = FALSE, dot.plot = TRUE)
  }
}

close.screen(all.screens=TRUE)
color.bar(colorRampPalette(c("blue", "grey", "red"))(100), 1)

plot(0, type='n', xlim=c(0, 5), ylim=c(0,2), axes=FALSE, ylab="", xlab="")
points(x = seq(0, 5), y=rep(1, 6), pch=16, 
       cex=c(rep(2, 2), seq(2,5)))
axis(side = 1, at=c(0:5), labels = c(0:5))

dev.off()

