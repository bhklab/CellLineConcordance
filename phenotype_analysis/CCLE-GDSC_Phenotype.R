
packages <- c("AnnotationDbi", "org.Hs.eg.db", "plyr",
             "vioplot", "PharmacoGx", "scales", "gtools",
             "intervals", "preprocessCore", "Hmisc", 
             "ggplot2", "parallel")
packages <- c("AnnotationDbi", "org.Hs.eg.db", "PharmacoGx", "scales")
tmp <- suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE))

outdir <- '~/Desktop/bhk_lab/results/phenotypes/expr_plots'
outdir <- file.path(outdir, "plots")
refdir <- '~/git/reference/l1000'
abc.dir <- '~/git/ccl_phenotype/reference'
src.dir <- '~/git/ccl_phenotype/src'

outdir <- '/results/pharmaco'  #codeocean
refdir <- '/data/ref'  #codeocean
abc.dir <- '/data/pharmaco'  #codeocean
src.dir <- '/code/src'  #codeocean

cloi <- "REH"

source(file.path(src.dir, "comparePhenotypes.R"))
load(file.path(refdir, "l1000.Rdata"))

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
  print(paste0("Tissue ID: ", getTissue(targ.cl, ds1)))
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
  if(mat.to.vector) matx <- as.numeric(round(matx,3))
  summ.metric <- list("Matching"=matching.diag,
                      "Nonmatching"=matx)
  if(mat.to.vector){
    if(return.style == 'data.frame' | return.style == 'melt'){
      require(rowr)
      dat <- cbind.fill(summ.metric[[1]], summ.metric[[2]])
      colnames(dat) <- c("Matching", "Nonmatching")
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
                            cex.range=c(2, 5), scale=1,
                            ramp.cols=list("1" = "black", "-1" = "red"),
                            ...){
  ks.pval <- ks.test(seq1, seq2)$p.val
  
  dseq1 <- density(seq1, na.rm=TRUE)
  dseq2 <- density(seq2, na.rm=TRUE)

  f1 <- approxfun(dseq1)
  f2 <- approxfun(dseq2)
    
  # Generates the KDE for the match/nonmatch curves
  if(kde.plot){   
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
      bf <- log2(p2/p1)
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
      
      p2 <- f2(x)
      p1 <- f1(x)
      if(!is.na(x) && is.na(p2)) p2 <- 0
      if(!is.na(x) && is.na(p1)) p1 <- 0
      bf <- round(log2(p2/p1), 1)
      
      if(!is.na(bf)){
        if(!is.list(ramp.cols)){
          ramp.cols <- list("1" = ramp.cols[1],
                            "-1" = ramp.cols[2])
        }
        cex.val <- abs(bf) * scale
        cex.range.r <- cex.range * scale
        if(cex.val > max(cex.range.r)) cex.val <- max(cex.range.r)
        if(cex.val < min(cex.range.r)) cex.val <- min(cex.range.r)
        
        points(x, each.col, pch=16, cex=cex.val, col=ramp.cols[[as.character(sign(bf))]])
      }
    }
  }

  return(list('ks'=ks.pval,
              'bf'=do.call("cbind", bf.vals)))
}

filterResistant <- function(AAClist, dsA, dsB, each.row, S.thr=0.2, Rtype='both', filter=TRUE){
    # Creates a 1xd matrix for gCSI drugs and dx1 matrix for CCLE drugs less than Sensitivity threshold
    row.R <- AAClist[[dsA]][each.row, ,drop=FALSE] < S.thr
    col.R <- AAClist[[dsB]][each.row, ,drop=FALSE] < S.thr

    # Remove resistant     if BOTH (2/2) are resistant
    R.mat <- t(row.R) %*% (col.R) # row.R  (n x 1) %*% col.R (1 x n)
    if(Rtype=='single'){
        # Remove resistant   if at least 1 of the 2 is resistant 
        R.mat[] <- 0L
        R.mat.single <- sweep(sweep(R.mat, 2, col.R, "+"), 1, t(row.R), "+") >= 1
        R.mat <- R.mat.single      
    }  
    mode(R.mat) <- "logical"
    if(filter){
        print(paste0("Removing: ",  rownames(AAClist[[dsA]])[each.row], " - ",
                     sum(R.mat, na.rm=TRUE), 
                     "/", sum(!is.na(R.mat))))
    }

    R.mat
}

computeDeltaABC <- function(AAClist, dsA, dsB, ABC, filter=TRUE, ...){
    Y.dABC <- lapply(seq_along(rownames(AAClist[[dsA]])), 
                               function(each.row, ...){                        
        # All-by-all comparison for ABC between gCSI (dsA) and CCLE (dsB)
        Y.tmp <- ABC[[rownames(AAClist[[dsA]])[each.row]]]
                                   
        # Creates a matrix of resistant lines
        R.mat <- filterResistant(AAClist, dsA, dsB, each.row, filter=filter, ...)
        c.idx <- match(colnames(R.mat), colnames(Y.tmp))
        r.idx <- match(rownames(R.mat), rownames(Y.tmp))
        Y.tmp <- Y.tmp[r.idx, c.idx]
        colnames(Y.tmp) <- rownames(Y.tmp) <- colnames(R.mat)
        if(!filter){
            R.mat <- matrix(NA)
        }


        # For all the Resistant lines (single or both), set that deltaAUC to NA
        diagY <- diag(Y.tmp)
        Y.tmp[t(R.mat)] <- NA
        diag(Y.tmp) <- diagY
        Y.tmp
    })

    names(Y.dABC) <- rownames(drug.auc.list[[dsA]])
    summ.abc <- lapply(Y.dABC, summMatchMatrix, mat.to.vector=TRUE, return.style='list')
    summ.abc              
}

computeDeltaAAC <- function(AAClist, dsA, dsB, filter=TRUE, ...){
    Y.dAAC <- lapply(seq_along(rownames(AAClist[[dsA]])), 
                               function(each.row, ...){                        
        # All-by-all comparison for delta-AUC between gCSI (dsA) and CCLE (dsB)
        Y.tmp <- apply(AAClist[[dsB]][each.row, ,drop=FALSE], 2, function(each.col){
            abs(each.col - AAClist[[dsA]][each.row, ])
        })
        # Assemble the matrix of deltaAUCs
        #Y.tmp <- sapply(Ylist.tmp, function(x) x[['dAUC']])
        rownames(Y.tmp) <- colnames(AAClist[[dsA]])
        colnames(Y.tmp) <- colnames(AAClist[[dsB]])

        # Creates a matrix of resistant lines
        if(filter){
            R.mat <- filterResistant(AAClist, dsA, dsB, each.row, ...)
        } else {
            R.mat <- NULL
        }


        # For all the Resistant lines (single or both), set that deltaAUC to NA
        diagY <- diag(Y.tmp)
        Y.tmp[R.mat] <- NA
        diag(Y.tmp) <- diagY
        Y.tmp
    })

    names(Y.dAAC) <- rownames(drug.auc.list[[dsA]])
    summ.aac <- lapply(Y.dAAC, summMatchMatrix, mat.to.vector=TRUE, return.style='list')
    summ.aac              
}

plotKs <- function(i, j, col.i='blue', col.j='green', col.D="red",
                   p=NA, D=NA,
                   include.d=TRUE, annotate=FALSE, add=FALSE, ...){
    # Form ecdf
    cdf1=ecdf(i)
    cdf2=ecdf(j)
    
    # Plot the ecdfs
    plot(cdf1, verticals=TRUE, do.points=FALSE, col=col.i, add=add, ...) 
    plot(cdf2, verticals=TRUE, do.points=FALSE, col=col.j, add=TRUE, ...) 
    
    # Mark the D
    if(include.d){
        minMax <- seq(min(i, j), max(i, j), length.out=length(i)) 
        x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )][1] 
        y0 <- cdf1(x0) 
        y1 <- cdf2(x0) 
        
        points(c(x0, x0), c(y0, y1), pch=16, col=col.D) 
        segments(x0, y0, x0, y1, col=col.D, lty="dotted") 
    }
    
    #Label the D and p val 
    if(annotate){
        if(!include.d){
            x0 <- min(c(i, j))
            y0 <- 1
            y1 <- 0.8
        }
        text(x=x0, y=(y0+y1)/2, labels = paste0("D=",D), adj = 0, ...)
        text(x=min(c(i,j)), y=1, labels = paste0("p=",p), adj=0,...)
    }
}

#### MAIN: Obtain PharmacoGX Data ####
#Download PSets
availablePSets()
GDSC <- downloadPSet("GDSC1000")
CCLE <-downloadPSet("CCLE")

#   Cell lines found in both GDSC and CCLE
commonCl <- intersectPSet(list('CCLE'=CCLE,
                               'GDSC'=GDSC),  
                          intersectOn=c("cell.lines"))

#   Cell lines found in both GDSC and CCLE and their same drug profiles
commonClAuc <- intersectPSet(list('CCLE'=CCLE,
                                  'GDSC'=GDSC), 
                             intersectOn=c("cell.lines", "drugs"))

#   Genes found in both GDSC and CCLE assay
commonGenes <- intersect(fNames(GDSC, "rna"),
                         fNames(CCLE,"rna"))

#   List of cell lines of the same tissue type as CLOI
tissueCls <- intersect(cellNames(commonClAuc[[1]]),
                       getTissueCommonCl(GDSC, CCLE, cloi))
commonTissueCl <- intersectPSet(list('CCLE'=CCLE,
                                     'GDSC'=GDSC),
                                intersectOn = c("drugs", "cell.lines"),
                                cells = tissueCls)

# EXPR for all cell lines in GDSC and CCLE
expr.list <- list()
expr.list[['CCLE.all']] <- getExprVals(CCLE, NULL, "rna", commonGenes)
expr.list[['GDSC.all']] <- getExprVals(GDSC, NULL, "rna", commonGenes)

# EXPR for cell lines found in both GDSC and CCLE
expr.list[['GDSC']] <- getExprVals(commonCl, 'GDSC', "rna", commonGenes)
expr.list[['CCLE']] <- getExprVals(commonCl, 'CCLE', "rna", commonGenes)

# EXPR for all cell lines of the same tissue type as CLOI
expr.list[['GDSC.tissue']] <- getExprVals(commonTissueCl, 'GDSC', "rna", commonGenes)
expr.list[['CCLE.tissue']] <- getExprVals(commonTissueCl, 'CCLE', "rna", commonGenes)

# Find the L1000 Genes Expression values for overlapping cell lines:
expr.sub.list <- list()
expr.sub.list[['GDSC']] <- subsetExprMat(exprs(expr.list[['GDSC']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)
expr.sub.list[['CCLE']] <- subsetExprMat(exprs(expr.list[['CCLE']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)

# L1000 gene expression for cell lines of same tissue type as CLOI
expr.sub.list[['GDSC.tissue']] <- subsetExprMat(exprs(expr.list[['GDSC.tissue']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)
expr.sub.list[['CCLE.tissue']] <- subsetExprMat(exprs(expr.list[['CCLE.tissue']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)

drug.auc.list <- list()
# Drug AUC for all cell lines in GDSc and CCLE
drug.auc.list[['CCLE.all']] <- getDrugVals(CCLE)
drug.auc.list[['GDSC.all']] <- getDrugVals(GDSC)

# Find the Drug Recomputed AUC values for overlapping cell lines:
drug.auc.list[['GDSC']] <- getDrugVals(commonClAuc, 'GDSC')
drug.auc.list[['CCLE']] <- getDrugVals(commonClAuc, 'CCLE')

# Drug AUC for cell lines of same tissue type as CLOI
drug.auc.list[['GDSC.tissue']] <- getDrugVals(commonTissueCl, 'GDSC')
drug.auc.list[['CCLE.tissue']] <- getDrugVals(commonTissueCl, 'CCLE')

#### Expression correlation matrices ####
# Pearson correlation of L1000 genes
Y.l1000 <- cor(as.matrix(expr.sub.list[['GDSC']]), 
               as.matrix(expr.sub.list[['CCLE']]), 
               use = "pairwise.complete.obs", method = "pearson")
summ.expr.l1000 <- summMatchMatrix(Y.l1000, TRUE, 'list')

# Pearson correlation of L1000 genes for tissue-specificy
Y.l1000 <- cor(as.matrix(expr.sub.list[['GDSC.tissue']]), 
               as.matrix(expr.sub.list[['CCLE.tissue']]), 
               use = "pairwise.complete.obs", method = "pearson")
summ.expr.l1000.tissue <- summMatchMatrix(Y.l1000, TRUE, 'list')

# All gene expression
Y <- cor(as.matrix(exprs(expr.list[['GDSC']])), 
         as.matrix(exprs(expr.list[['CCLE']])), 
         use = "pairwise.complete.obs", method = "pearson")
summ.expr <- summMatchMatrix(Y, TRUE, 'list')
 
# All gene expression in same tissue type as CLOI
Y <- cor(as.matrix(exprs(expr.list[['GDSC.tissue']])), 
         as.matrix(exprs(expr.list[['CCLE.tissue']])), 
         use = "pairwise.complete.obs", method = "pearson")
summ.expr.tissue  <- summMatchMatrix(Y, TRUE, 'list')

#### Drug sensitivity matrices ####
load(file.path(abc.dir, "abc.gdsc-ccle.Rdata"))

#abc.matrices[sapply(abc.matrices, class) == 'try-error'] <- NULL
dABC.all <- computeDeltaABC(drug.auc.list, dsA='GDSC', dsB='CCLE', ABC=abc.matrices,
                            filter=FALSE)
S.dABC.all <- computeDeltaABC(drug.auc.list, dsA='GDSC', dsB='CCLE', ABC=abc.matrices,
                              filter=TRUE, S.thr=0.2, Rtype='both')
paclitaxel <- computeDeltaABC(drug.auc.list, dsA='GDSC', dsB='CCLE',  ABC=abc.matrices,
                              filter=TRUE, S.thr =0.4, Rtype='both')[['paclitaxel']]
S.dABC.all[['paclitaxel']] <- paclitaxel


#### Visualization: CDFs for Expression ####
out.id <- 'paper'

switch(out.id,
       reh={
         ids.x <- 'REH'
       },
       paper={
         ids.x <- c("REH", "NCI-H23")  
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


pdf(file.path(outdir, "ks_l1000Expr.pdf"), height=3, width=15)
close.screen(all.screens=TRUE)
split.screen(c(1,length(S.dABC.all)))
I <- summ.expr[['Matching']]
J <- summ.expr[['Nonmatching']]
KS <- suppressWarnings(ks.test(I, J))

screen(1)
par(mar=c(5.1,0.5,4.1,0.5))
p <- round(KS$p.value,4)
D <- round(KS$statistic, 3)

plot(0, type='n', xlim=c(0,1), ylim=c(0,1), main="Expression", 
     xlab=paste(c(D,if(p > 0) p else "<0.0001"), collapse="\n"),
     axes=FALSE, ylab='', sub='', cex.main=0.8, cex.sub=0.8)
plotKs(na.omit(I), na.omit(J), add=TRUE, include.d = TRUE,
       col.i=alpha("black", 0.8), col.j=alpha("red", 0.8), col.D="black",
       p=p, D=D, annotate=FALSE)
axis(side=1, at = c(0,1), labels = c("0", "1"))

dev.off()

#### Visualization: CDFs for Drug sensitivity ABC ####
print("KS Test for Matching/Nonmatching drug ABC curves")
print(sapply(names(dABC.all), function(drug.id){
    round(suppressWarnings(
      ks.test(dABC.all[[drug.id]][['Matching']], 
              dABC.all[[drug.id]][['Nonmatching']])$p.val), 5)
}))



## Resistant removed ABCs
pdf(file.path(outdir, "ks_drugABC.pdf"), height=3, width=15)
close.screen(all.screens=TRUE)
split.screen(c(1,length(S.dABC.all)))
kspval <- sapply(names(S.dABC.all), function(drug.id, plot.ks){
    I <- S.dABC.all[[drug.id]][['Matching']]
    J <- S.dABC.all[[drug.id]][['Nonmatching']]
    KS <- suppressWarnings(ks.test(I, J))

    if(plot.ks){
        screen(grep(drug.id, names(S.dABC.all)))
        par(mar=c(5.1,0.5,4.1,0.5))
        p <- round(KS$p.value,4)
        D <- round(KS$statistic, 3)
        
        plot(0, type='n', xlim=c(0,1), ylim=c(0,1), main=drug.id, 
             xlab=paste(c(D,if(p > 0) p else "<0.0001"), collapse="\n"),
             axes=FALSE, ylab='', sub='', cex.main=0.8, cex.sub=0.8)
        plotKs(na.omit(I), na.omit(J), add=TRUE, include.d = TRUE,
               col.i=alpha("black", 0.8), col.j=alpha("red", 0.8), col.D="black",
               p=p, D=D, annotate=FALSE)
        axis(side=1, at = c(0,1), labels = c("0", "1"))
    }
    
    pval <- KS$p.val
    

    return(c("meanM"=mean(I, na.rm=TRUE),
        "sdM"=sd(I, na.rm=TRUE),
        "meanNM"=mean(J, na.rm=TRUE),
        "sdNM"=sd(J, na.rm=TRUE),
        "pval"=pval))
}, plot.ks=TRUE)
dev.off()


meanM <- t(kspval[c(1,3),])
colnames(meanM) <- c("Matching", "Non-matching")

pdf(file.path(outdir, "avg_deltaABC.pdf"))
box.points <- boxplot(meanM, las=1,
                      ylab="Average Delta ABC",
                      ylim=c(0, 0.5), main="Resistant-removed Delta ABC")
points(rep(1, nrow(meanM)), meanM[,1], col="black", pch=16)
text(rep(1.1, nrow(meanM)), meanM[,1], labels=rownames(meanM), adj = 0, cex=0.7)
points(rep(2, nrow(meanM)), meanM[,2], col="black", pch=16)
text(rep(2.1, nrow(meanM)), meanM[,2], labels=rownames(meanM), adj = 0, cex=0.7)
dev.off()



plot.pdf <- TRUE

if(plot.pdf) pdf(file.path(outdir, "density_abc.reh-nci.pdf"), height=5, width=5)
lr.abc <- sapply(names(S.dABC.all), function(drug.id){
    values.x <- getXVals(S.dABC.all[[drug.id]][[1]], ids.x)

    boxplot(S.dABC.all[[drug.id]])
    plotBayesFactor(S.dABC.all[[drug.id]][[1]], S.dABC.all[[drug.id]][[2]],
                  values.x = values.x, f1col = "black",
                  xlab='ABC', main=drug.id)
})
lr.abc.df <- do.call("rbind", lr.abc['bf',])
rownames(lr.abc.df) <- colnames(lr.abc)
lr.abc.df
if(plot.pdf) dev.off()

