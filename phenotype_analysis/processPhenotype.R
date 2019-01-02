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

source("~/git/ccl_phenotype/src/comparePhenotypes.R")
load("/Users/rquevedo/git/reference/l1000/l1000.Rdata")
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

#### Obtain all PharmacoGX Data
#Download PSets
availablePSets()
GDSC <- downloadPSet("GDSC")
load("~/Desktop/bhk_lab/data/pharmacogx/GDSC.RData")
CCLE <- downloadPSet("CCLE")

commonClAuc <- intersectPSet(list('CCLE'=CCLE,
                             'GDSC'=GDSC), 
                             intersectOn=c("cell.lines", "drugs"))
commonCl <- intersectPSet(list('CCLE'=CCLE,
                             'GDSC'=GDSC), 
                          intersectOn=c("cell.lines"))
commonGenes <- intersect(fNames(GDSC, "rna"),
                         fNames(CCLE,"rna"))

# Find the Expression values for overlapping cell lines:
expr.list <- list()
expr.list[['CGP']] <- summarizeMolecularProfiles(commonCl$GDSC,
                                             cellNames(commonCl$GDSC),
                                             mDataType="rna",
                                             features=commonGenes,
                                             verbose=FALSE)
expr.list[['CCLE']] <- summarizeMolecularProfiles(commonCl$CCLE,
                                             cellNames(commonCl$CCLE),
                                             mDataType="rna",
                                             features=commonGenes,
                                             verbose=FALSE)
expr.list[['CCLE.all']] <- summarizeMolecularProfiles(CCLE,
                                                  cellNames(CCLE),
                                                  mDataType="rna",
                                                  features=commonGenes,
                                                  verbose=FALSE)
expr.list[['GDSC.all']] <- summarizeMolecularProfiles(GDSC,
                                                      cellNames(GDSC),
                                                      mDataType="rna",
                                                      features=commonGenes,
                                                      verbose=FALSE)
expr.list[['GDSC']] <- summarizeMolecularProfiles(commonCl$GDSC,
                                              cellNames(commonCl$GDSC),
                                              mDataType="rna2",
                                              features=commonGenes,
                                              verbose=FALSE)
#x <- cbind(exprs(expr.list[['CCLE.all']])[,'OCI-LY10'], exprs(expr.list[['GDSC.all']])[,'KARPAS-422'])
#x <- cbind(exprs(expr.list[['CCLE.all']])[,'SK-BR-3'], exprs(expr.list[['CCLE.all']])[,'AU565'])

# Subset the Expr data for the L1000 Genes
expr.sub.list <- list()
expr.sub.list[['CGP']] <- subsetExprMat(exprs(expr.list[['CGP']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)
expr.sub.list[['GDSC']] <- subsetExprMat(exprs(expr.list[['GDSC']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)
expr.sub.list[['CCLE']] <- subsetExprMat(exprs(expr.list[['CCLE']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)


# Find the Drug Recomputed AUC values for overlapping cell lines:
drug.auc.list <- list()
drug.auc.list[['GDSC']] <- summarizeSensitivityProfiles(
                              pSet=commonClAuc$GDSC,
                              sensitivity.measure='auc_recomputed',
                              summary.stat='median')
drug.auc.list[['CCLE']] <- summarizeSensitivityProfiles(
                              pSet=commonClAuc$CCLE,
                              sensitivity.measure='auc_recomputed',
                              summary.stat='median')
drug.auc.list[['CCLE.all']] <- summarizeSensitivityProfiles(
                              pSet=CCLE,
                              sensitivity.measure='auc_recomputed',
                              summary.stat='median')
drug.auc.list[['GDSC.all']] <- summarizeSensitivityProfiles(
                              pSet=GDSC,
                              sensitivity.measure='auc_recomputed',
                              summary.stat='median')

#### Process the mismatch, match, and concordance IDs
getCcleGdscCnvDiff <- function(x){
  ccle.row <- grep("ccle", rownames(x))
  ccle.cnvDiff <- x[ccle.row, (x[ccle.row,] < 0.8), drop=FALSE]
  col.idx <- sapply(c("cgp", "gdsc"), function(x) grep(x, colnames(ccle.cnvDiff)))
  col.idx <- unlist(col.idx)
  ccle.cnvDiff <- ccle.cnvDiff[,col.idx, drop=FALSE]
  return(ccle.cnvDiff)
}

# Cell lines with same annotations but discordant SNP genotypes
load("/Users/rquevedo/git/snp_fingerprint/environments/ccle.matchAnno.Rdata") # ccle.matchAnno.df
snpDiff.cellids <- which(apply(ccle.matchAnno.df, 1, function(x) any(x['cgp'] < 1 | x['gdsc'] < 1)))

# Cell lines with same annotations but discordant CNV profiles
load("/Users/rquevedo/git/snp_fingerprint/environments/cnv.disc.match.Rdata") # cnv.disc.match
cnvDiff.list <- list()
cnvDiff.list[['hscra1']] <- lapply(cnv.disc.match[['hscra1']][['ccle']], getCcleGdscCnvDiff)
cnvDiff.list[['hscra2']] <- lapply(cnv.disc.match[['hscra2']][['ccle']], getCcleGdscCnvDiff)

# Cell lines with different annotations but concordant genotypes
# [1] "matrix.perc.match"            "mm.ds.list"                   "unknown.ds.list"              "mismatch.m.filt.df"
# [5] "mismatch.m.df.bkup"           "combined.match.mismatch.list"
load("/Users/rquevedo/git/snp_fingerprint/environments/subHeatmap.Feb17-2016.Rdata")
formatMismatch <- function(x){
  ds.ids <- sapply(x, function(y) gsub("_.+", "", y))
  cl.ids <- sapply(x, function(y) gsub("^.+?_", "", y, perl=TRUE))
  
  ds.cl.idx <- which(ds.ids %in% c('ccle', 'gdsc', 'cgp'))
  ds.ids <- ds.ids[ds.cl.idx]
  cl.ids <- cl.ids[ds.cl.idx]
  return(list("ds"=ds.ids,
              "cl"=cl.ids))
}

snp.mismatch.df <- do.call("rbind.fill", mm.ds.list[c("ccle", "gdsc.x", "cgp")])
snp.mismatch.df <- as.data.frame(gsub("gdsc.x", "gdsc", as.matrix(snp.mismatch.df)))
snp.mismatch.list <- apply(snp.mismatch.df, 1, formatMismatch)


#### -----------------------------------------------------------------
#### Compute the ECDFs for Deltas in drug AUC between cell lines
getDeltaAuc <- function(auc.x, auc.y, cell.ids, 
                        analysis.type='concordant', filt=FALSE,
                        sub.size=1.0, targ.cl='any'){
  if(nrow(auc.x) == nrow(auc.y)){ # Number of drugs is the same between X and Y
    if(analysis.type == 'concordant'){  # Cell lines should be the same
      delta.auc.list <- list()
      sub.cell.ids <- sample(x = cell.ids, size = (sub.size * length(cell.ids)))
      while(!targ.cl %in% sub.cell.ids){
        sub.cell.ids <- sample(x = cell.ids, 
                               size = (sub.size * length(cell.ids)))
        if(targ.cl=='any') targ.cl <- sub.cell.ids[1]
      }
      auc.x <- auc.x[,sub.cell.ids, drop=FALSE]
      auc.y <- auc.y[,sub.cell.ids, drop=FALSE]
      
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

cnv.ids <- unique(sort(unlist(lapply(cnvDiff.list, names))))
cnv.ids <- c("HT-29", "HuP-T4", "MCF7", "Mewo", "NCI-H23", "REH", "UACC-812")
cnv.ids <- c("REH")
cnv.ids <- c("KNS-81-FD")
# Delta AUC for matching cell lines that have discordant CNVs
cnvDiff.delta.auc <- getDeltaAuc(auc.x=drug.auc.list[['CCLE']], 
                                 auc.y=drug.auc.list[['GDSC']], 
                                 cell.ids=cnv.ids, analysis.type='concordant',
                                 filt=FALSE)
# Delta AUC for matching cell lines that have discordant SNPs
snp.ids <- names(snpDiff.cellids)
snpDiff.delta.auc <- getDeltaAuc(auc.x=drug.auc.list[['CCLE']], 
                                 auc.y=drug.auc.list[['GDSC']], 
                                 cell.ids=snp.ids, analysis.type='concordant',
                                 filt=FALSE)
all.cell.ids <- c(cnv.ids, snp.ids)
ref.ids <- colnames(drug.auc.list[['CCLE']])[-which(colnames(drug.auc.list[['CCLE']]) %in% all.cell.ids)]
# Delta AUC for all matching cell lines
delta.auc <- getDeltaAuc(auc.x=drug.auc.list[['CCLE']], 
                         auc.y=drug.auc.list[['GDSC']], 
                         cell.ids=ref.ids, analysis.type='concordant',
                         filt=TRUE)
# Delta AUC for all nonmatching cell lines
disc.delta.auc <- getDeltaAuc(auc.x=drug.auc.list[['CCLE']], 
                              auc.y=drug.auc.list[['GDSC']], 
                              cell.ids=NA, analysis.type='discordant',
                              filt=TRUE)

#### -----------------------------------------------------------------
#### ------------------Permutation Testing----------------------------
deltaAucPermutation <- function(targ.cl='any', nperm=1000, sub.size=0.6, 
                                ref.ids, auc.list, top.val=0.05,
                                analysis.type='concordant', 
                                ds1='CCLE', ds2='GDSC',
                                cl1=NULL, cl2=NULL, 
                                alt.ds1=NULL, alt.ds2=NULL){
  ## Wrapper for the **getDeltaAuc()** function that will do n number of
  ## permutations and find out how many times a target cell line is found
  ## in the top X of all discordant cell lines for that drug
  targ.cl.tmp <- targ.cl
  top.val.cnt <- sapply(seq(nperm), function(x){
    # Get the Delta AUC for each cell line for each drug
    delta.auc <- getDeltaAuc(auc.x=auc.list[[ds1]], auc.y=auc.list[[ds2]], 
                             cell.ids=ref.ids, analysis.type=analysis.type,
                             filt=FALSE, sub.size=0.6, targ.cl=targ.cl)
    
    # Checks to see if two cell lines with different anno's were given and sets those
    # as the target adter appending
    if(!is.null(cl1) && !is.null(cl2)){
      # Create a dataframe containing just the two cell lines of interest for overlapping drugs
      alt.ds2 <- 'CCLE.all'
      auc.x <- auc.list[[alt.ds1]]
      auc.y <- auc.list[[alt.ds2]]
      auc.x <- auc.x[which(rownames(auc.x) %in% rownames(delta.auc)),]
      auc.y <- auc.y[which(rownames(auc.y) %in% rownames(delta.auc)),]

      auc.xy = data.frame(cl1=auc.x[,cl1, drop=FALSE],
                          cl2=auc.y[,cl2, drop=FALSE])
      delta.auc.xy <- as.matrix(abs(apply(auc.xy, 1, diff)), ncol=1)
      colnames(delta.auc.xy) <- paste(cl1, cl2, sep="_")
      
      if(all(rownames(delta.auc.xy) == rownames(delta.auc))){
        delta.auc <- cbind(delta.auc, delta.auc.xy)
        targ.cl.tmp <- colnames(delta.auc.xy)
      } else {
        stop("Rownames were not ordered the same!")
      }
    }
    
    # Order the delta-auc's so the most discordant for each drug is the first in the list
    apply(delta.auc, 1, function(drug.x.auc, targ.cl){
      ord.auc <- colnames(delta.auc)[rev(order(drug.x.auc,na.last=FALSE))]
      #TRUE if its within the top X most discordant AUCs
      val.cutoff <- ceiling(length(na.omit(drug.x.auc)) * top.val)
      if(val.cutoff < val.cutoff) val.cutoff <- 3
      grep(paste0("^", targ.cl, "$"), ord.auc) <= val.cutoff
    }, targ.cl=targ.cl.tmp)
    
  })
  
  print("Getting pval")
  pval <- as.matrix((1 - (apply(top.val.cnt, 1, sum) / nperm)), ncol=1)
  colnames(pval) <- targ.cl.tmp
  pval
}
reh.pval <- deltaAucPermutation(targ.cl='REH', nperm=1000, sub.size=0.6,
                                top.val=0.1,
                                ref.ids=colnames(drug.auc.list[['CCLE']]), 
                                auc.list=drug.auc.list,
                                ds1='CCLE', ds2='GDSC')
nci.pval <- deltaAucPermutation(targ.cl='NCI-H23', nperm=1000, sub.size=0.6,
                                top.val=0.1,
                                ref.ids=colnames(drug.auc.list[['CCLE']]), 
                                auc.list=drug.auc.list,
                                ds1='CCLE', ds2='GDSC')
oci_karpas.pval <- deltaAucPermutation(targ.cl='any', nperm=1000, 
                                       sub.size=0.6, top.val=0.1,
                                ref.ids=colnames(drug.auc.list[['CCLE']]), 
                                auc.list=drug.auc.list,
                                ds1='CCLE', ds2='GDSC',
                                cl1='OCI-LY10', cl2='KARPAS-422',
                                alt.ds1='CCLE.all', alt.ds2='CCLE.all')
kns81fd.pval  <- deltaAucPermutation(targ.cl='KNS-81-FD', top.val=0.1,
                                     nperm=100, sub.size=0.6,
                                ref.ids=colnames(drug.auc.list[['CCLE']]), 
                                auc.list=drug.auc.list,
                                ds1='CCLE', ds2='GDSC')



pdf("/Users/rquevedo/Desktop/bhk_lab/results/phenotypes/summary_metrics/drug_auc.pdf")
scr.spl <- ceiling(sqrt(nrow(delta.auc)))
split.screen(c(scr.spl, scr.spl))
for(i in c(1:nrow(delta.auc))){
  screen(i)
  par(mar=c(2, 2, 2, 2))
  drug.id <- rownames(delta.auc)[i]
  tryCatch(plot(ecdf(delta.auc[drug.id,]), main = drug.id, 
                do.points=FALSE, verticals=TRUE, xlab="", xlim=c(0,1), las=2, cex.axis=0.7),
           error=function(e){plot(0, type='n', xlim=c(0,1), ylim=c(0,1))})
  disc.drug.id.idx <- match(drug.id, rownames(disc.delta.auc[[2]]))
  disc.drug.ecdf <- sapply(disc.delta.auc, function(x) {
    tryCatch(as.numeric(x[disc.drug.id.idx,]),
             error=function(e){NA})
    })
  tryCatch(plot(ecdf(disc.drug.ecdf), 
                do.points=FALSE, verticals=TRUE, add=TRUE, col="grey"),
           error=function(e){print(paste("No SNP Discordance for", drug.id))})
  tryCatch(plot(ecdf(snpDiff.delta.auc[drug.id,]), 
                do.points=FALSE, verticals=TRUE, add=TRUE, col="magenta"),
           error=function(e){print(paste("No SNP Discordance for", drug.id))})
  tryCatch(plot(ecdf(cnvDiff.delta.auc[drug.id,]), 
                do.points=FALSE, verticals=TRUE, add=TRUE, col="blue"),
           error=function(e){print(paste("No CNV Discordance for", drug.id))})
}
screen(16)
par(mar=c(2, 2, 2, 2))
plot(0, type='n', xlim=c(0,10), ylim=c(1,9), axes=FALSE, xlab='', ylab='')
rect(xleft = rep(1, 4), ybottom = seq(2,8, by=2), 
     xright = rep(3, 4), ytop = seq(3,9, by=2),
     col = c("black", "grey", "magenta", "blue"))
text(x = rep(3,4), y=seq(2.5, 8.5, by=2), 
     labels=c("Ref match", "Ref disc", "snp disc", "cnv disc"), pos=4, cex=0.5)
close.screen(all.screens=TRUE)
dev.off()



#### Map the Delta AUCs to the reference ECDFS between cell lines
getCnvDisc <- function(each.hscr, cnv.ids){
  unlist(sapply(cnv.ids, function(y) {
    x <- cnvDiff.list[[each.hscr]][[y]]
    if(is.null(x)){ return(NA) }
    if(ncol(x) == 0){ return(NA) }
    row.idx <- grep("gdsc", colnames(x))
    if(length(row.idx) == 0) row.idx <- 1
    return(x[,row.idx])
  }))
}


pdf("/Users/rquevedo/Desktop/bhk_lab/results/phenotypes/summary_metrics/deltaAuc_density.kns-81-fd.pdf")
split.screen(c(3,2))
screen.id <- 1
drug.list <- c("AZD6244", "PD-0332991", "PD-0325901", "paclitaxel", "17-AAG")
drug.list <- c("PD-0332991", "paclitaxel", "17-AAG")
for(drug.id in drug.list){  
  screen(screen.id)
  par(mar=c(2,2,2,2))
  # For "drug.id", pulls out the delta AUC for every possible cell line pairing 
  disc.drug.id.idx <- match(drug.id, rownames(disc.delta.auc[[2]]))
  disc.drug.ecdf <- sapply(disc.delta.auc, function(x) {
    tryCatch(as.numeric(x[disc.drug.id.idx,]),
             error=function(e){NA})
  })
  
  
  conc.dens.ecdf <- ecdf(delta.auc[drug.id,])  # Delta AUC for matching cell lines
  disc.dens.ecdf <- ecdf(disc.drug.ecdf)  # Delta AUC for all paired cell lines
  dens.conc.auc <- density(delta.auc[drug.id,], na.rm=TRUE)
  plot(dens.conc.auc, ylim=c(0,10), xlim=c(0,1), xlab="Delta AUC", main=drug.id)
  name.ids <- names(cnvDiff.delta.auc[drug.id,,drop=FALSE])
  if(is.null(name.ids)) name.ids <- colnames(cnvDiff.delta.auc[drug.id,,drop=FALSE])
  sapply(name.ids, function(auc.name){
    each.auc <- cnvDiff.delta.auc[drug.id,auc.name]
    xval <- round(each.auc, 2)
    yval <- mean(with(dens.conc.auc, y[which(round(x, 2) == round(each.auc, 2))]))
    if(is.finite(xval) & !is.finite(yval)) yval <- 0
    
    points(xval, yval)
    text(xval, (yval + 0.5), labels=auc.name, cex=0.8, pos=4, srt=70)
    })
  
  
  lines(density(disc.drug.ecdf, na.rm=TRUE), col="red")
  perc.cnv.ecdf <- sapply(cnvDiff.delta.auc[drug.id,,drop=FALSE], function(x) conc.dens.ecdf(x))
  perc.cnv.disc.ecdf <- sapply(cnvDiff.delta.auc[drug.id,,drop=FALSE], function(x) disc.dens.ecdf(x))
  perc.cnv.ecdf.df <- rbind(perc.cnv.ecdf, perc.cnv.disc.ecdf, cnvDiff.delta.auc[drug.id,,drop=FALSE], 
                            getCnvDisc("hscra1", cnv.ids), getCnvDisc("hscra2", cnv.ids))
  rownames(perc.cnv.ecdf.df) <- c("matching_ecdf_perc", "nonmatch_ecdf_perc", "delta_AUC", "hscra1_cnv_corr", "hscra2_cnv_corr")
  print(perc.cnv.ecdf.df)
  # perc.cnv.ecdf.df for cell line X:
  #   Percentile that the delta AUC for cell line X falls into for all matching cell-line deltaAUC
  #   Percentile that the delta AUC for cell line X falls into for all non-matching cell-lines deltaAUC
  #   Absolute value for cell line X delta AUCs
  #   CNV concordance for cell line X HSCRA1
  #   CNV concordance for cell line X HSCRA2
  screen.id <- screen.id + 1
}
close.screen(all.screens=TRUE)
dev.off()

#### Compute the concordance between Gene Expression profiles (L1000)
getCorrMatrix <- function(expr.x, expr.y, cell.ids, analysis.type='concordant'){
  if(nrow(expr.x) == nrow(expr.y)){
    if(analysis.type=="concordant"){
      # Generate correlation of all expressions between Dataset X and Dataset Y
      # for cell lines that are annotated the same 
      print("Obtaining gene expression correlation matrix for similarly annotated cell lines")
      expr.x <- expr.x[,intersect(cell.ids, colnames(expr.x))]
      expr.y <- expr.y[,intersect(cell.ids, colnames(expr.y))]
      expr.ds <- sapply(colnames(expr.x), function(each.cl){
        combine.expr <- cor(expr.x[, each.cl, drop=FALSE], 
                            expr.y[, each.cl, drop=FALSE],
                            method="pearson", use="pairwise.complete.obs")
        return(combine.expr)
      })
    } else if(analysis.type=="discordant"){
      # Generate correlation of all expressions between Dataset X and Dataset Y
      # for cell lines that are annotated differently. 
      print("Obtaining gene expression correlation matrix for differently annotated cell lines")
      if(is.na(cell.ids)){
        expr.ds <- sapply(colnames(expr.x), function(each.cl.x){
          combine.expr.row <- sapply(colnames(expr.y), function(each.cl.y){
            if(each.cl.x != each.cl.y){
              combine.expr <- cor(expr.x[, each.cl.x, drop=FALSE], 
                                  expr.y[, each.cl.y, drop=FALSE],
                                  method="pearson", use="pairwise.complete.obs")
            } else {
              combine.expr <- NA
              print(paste(each.cl.x, "matches", each.cl.y))
            }
            return(combine.expr)
          })
          return(combine.expr.row)
        })
        
      } else {
        expr.ds <- lapply(cell.ids, function(each.mm){
          expr.ref <- NA
          print(paste("---", each.mm[['ds']]['Ref'], each.mm[['cl']]['Ref'], sep="_"))
          tryCatch(if(each.mm[['ds']]['Ref'] %in% c("ccle")){
              expr.ref <- expr.x[,each.mm[['cl']]['Ref'], drop=FALSE]
            } else {
              expr.ref <- expr.y[,each.mm[['cl']]['Ref'], drop=FALSE]
            },
            error=function(e){print(paste("Could not find ref:", each.mm[['cl']]['Ref']))})
          
          combine.expr <- NA
          if((length(each.mm[['ds']]) > 1) & !is.na(expr.ref)){
            
            combine.expr <- lapply(c(2:length(each.mm[['ds']])), function(each.cl){
              expr.alt <- NA
              tryCatch(if(each.mm[['ds']][each.cl] %in% c("ccle")) {
                expr.alt <- expr.x[,each.mm[['cl']][each.cl], drop=FALSE]
              } else {
                expr.alt <- expr.y[,each.mm[['cl']][each.cl], drop=FALSE]
              },
              error=function(e){print(paste("Could not find alt:", each.mm[['cl']][[each.cl]]))})
              
              if(!is.na(expr.alt) & !is.na(expr.ref)){
                combine.expr <- cor(expr.ref, expr.alt,
                                    method="pearson", use="pairwise.complete.obs")
              }
            
              return(combine.expr)
            })
            combine.expr <- do.call("cbind", combine.expr)
          }
          
          return(combine.expr)
        })
      }
    } else {
      print("Please specify analysis.type as either 'concordant' or 'discordant'")
    }
  } else {
    print("Error: Intersect based on common cell lines first")
    expr.ds <- NA
  }
  return(expr.ds)
}

plotExprDensity <- function(expr.x.df, expr.y.df, snp.ids, cnv.ids, ref.ids,
                            cnv.mismatch.list, snp.mismatch.list){
  expr.ds <- getCorrMatrix(expr.x=expr.x.df, expr.y=expr.y.df, 
                           cell.ids=ref.ids, analysis.type='concordant')
  cnv.ds <- getCorrMatrix(expr.x=expr.x.df, expr.y=expr.y.df, 
                          cell.ids=cnv.ids, analysis.type='concordant')
  snp.ds <- getCorrMatrix(expr.x=expr.x.df, expr.y=expr.y.df, 
                          cell.ids=snp.ids, analysis.type='concordant')
  cnv.snp.ds <- getCorrMatrix(expr.x=expr.x.df, expr.y=expr.y.df, 
                              cell.ids=intersect(snp.ids, cnv.ids), 
                              analysis.type='concordant')
  expr.disc.ds <- getCorrMatrix(expr.x=expr.x.df, expr.y=expr.y.df, 
                                cell.ids=NA, analysis.type='discordant')
  # cnv.disc.ds <- getCorrMatrix(expr.x=expr.x.df, expr.y=expr.y.df, 
  #                              cell.ids=cnv.mismatch.list, analysis.type='discordant')
  snp.disc.ds <- getCorrMatrix(expr.x=expr.x.df, expr.y=expr.y.df, 
                               cell.ids=snp.mismatch.list, analysis.type='discordant')
 
  rescaleDens <- function(dens.func){
    dens.y <- dens.func$y
    dens.y <- rescale(dens.y, to=c(0,1))
    dens.func$y <- dens.y
    return(dens.func)
  }
  
  # x=1: Correlation for all concordant and didscordant annotated cell lines across the datasets
  vioplot(na.omit(expr.ds), c(-2,-1), c(-2,-1), ylim=c(0,1), col=alpha("grey", 0.3))
  vioplot(na.omit(as.vector(expr.disc.ds)), at=1, add=TRUE, col=alpha("grey", 0.3))
  text(x=1, y = 1, labels=paste("Concordant Anno: n=", length(na.omit(expr.ds)), sep=""), cex=0.7)
  text(x=1, y = 0.9, labels=paste("Discordant Anno: n=", length(na.omit(as.vector(expr.disc.ds))), sep=""), cex=0.7)
  # x=2: Correlation for all CNV discordant, similar annotation (blue) and SNP discordant, similar annotation (pink)
  vioplot(na.omit(cnv.ds), col=alpha("blue", 0.3), ylim=c(0,1), at=2, add=TRUE)
  vioplot(na.omit(snp.ds), col=alpha("magenta", 0.3), at=2, add=TRUE)
  text(x=2, y = 1, labels=paste("SNP disc (magenta): n=", length(na.omit(snp.ds)), sep=""), cex=0.7)
  text(x=2, y = 0.9, labels=paste("CNV disc (blue): n=", length(na.omit(cnv.ds)), sep=""), cex=0.7)
  # x=3: Correlation for all cell lines with different annotations but concordant genotype
  vioplot(na.omit(unlist(snp.disc.ds)), col=alpha("magenta", 0.3), at=3, add=TRUE)
  text(x=3, y = 1, labels=paste("SNP conc: n=", length(na.omit(unlist(snp.disc.ds))), sep=""), cex=0.7)
  
  # plot(rescaleDens(density(expr.ds, na.rm=TRUE)),ylim=c(0,1), xlim=c(0, 1))
  # lines(rescaleDens(density(expr.disc.ds, na.rm=TRUE)), col="black")
  # lines(rescaleDens(density(cnv.ds, na.rm=TRUE)), col="blue")
  # lines(rescaleDens(density(snp.ds, na.rm=TRUE)), col="magenta")
  # 
  # plot(rescaleDens(density(expr.disc.ds, na.rm=TRUE)), ,ylim=c(0,1), xlim=c(0, 1), col="black")
  # lines(rescaleDens(density(expr.ds, na.rm=TRUE)), col="black")
  # lines(rescaleDens(density(unlist(snp.disc.ds), na.rm=TRUE)), col="magenta")
  
  #lines(rescaleDens(density(cnv.snp.ds, na.rm=TRUE)), col="green")
  # COLO-783                    COR-L51 Ishikawa (Heraklio) 02 ER- 
  #   NA                         NA                         NA 
  # MOG-G-CCM                        NB4                      SW403 
  #   NA                         NA                         NA 
}

pdf("/Users/rquevedo/Desktop/bhk_lab/results/phenotypes/summary_metrics/exprDensPlots.pdf")
plotExprDensity(expr.x.df=expr.sub.list[['CCLE']],
                expr.y.df=expr.sub.list[['CGP']],
                snp.ids, cnv.ids, ref.ids,
                NA, snp.mismatch.list)
plotExprDensity(expr.x.df=exprs(expr.list[['CCLE']]),
                expr.y.df=exprs(expr.list[['CGP']]),
                snp.ids, cnv.ids, ref.ids,
                NA, snp.mismatch.list)
dev.off()


#### Generate gene-level biological variation (internal noise) models
validateGeneOrder <- function(expr.x.df, expr.y.df){
  # Check to make sure all rownames in x are found in y
  if(!all(rownames(expr.x.df) == rownames(expr.y.df))){
    rows.x <- rownames(expr.x.df)
    rows.y <- rownames(expr.y.df)
    
    intersect.rows <- intersect(rows.x, rows.y)
    expr.x.df <- expr.x.df[match(intersect.rows, rows.x),]
    expr.y.df <- expr.y.df[match(intersect.rows, rows.y),]
  }
  return(list("x"=expr.x.df, "y"=expr.y.df))
}


generateGeneDist <- function(expr.xy.list, type='match'){
  row.idx <- c(1:nrow(expr.xy.list[['x']]))
  cell.ids <- colnames(expr.xy.list[['x']])
  
  if(type == 'match'){
    match.drug.cl.vals <- lapply(row.idx, function(each.gene.idx) {
      match.cl.vals <- sapply(cell.ids, function(each.cl){
        x.gene.cl.val <- expr.xy.list[['x']][each.gene.idx, each.cl]
        y.gene.cl.val <- expr.xy.list[['y']][each.gene.idx, each.cl]
        diff.xy.val <- abs(x.gene.cl.val - y.gene.cl.val)
        return(diff.xy.val)
      })
      return(match.cl.vals)
    })
    comp.cl.vals <- match.drug.cl.vals
    
  } else if (type == 'mismatch'){
    mismatch.drug.cl.vals <- lapply(row.idx, function(each.gene.idx) {
      mismatch.cl.vals <- sapply(cell.ids, function(each.cl.x){
        diff.xy.vals <- sapply(cell.ids, function(each.cl.y){
          if(each.cl.x != each.cl.y){
            x.gene.cl.val <- expr.xy.list[['x']][each.gene.idx, each.cl.x]
            y.gene.cl.val <- expr.xy.list[['y']][each.gene.idx, each.cl.y]
            diff.xy.val <- abs(x.gene.cl.val - y.gene.cl.val)
            return(diff.xy.val)
          }
        })
        return(unlist(diff.xy.vals))
      })
    })

    comp.cl.vals <- mismatch.drug.cl.vals
    
  } else {
    print("Arguement 'type' must be either 'match' or 'mismatch'")
    comp.cl.vals <- NA
  }
  
  names(comp.cl.vals) <- rownames(expr.xy.list[['x']])
  return(comp.cl.vals)
}

expr.xy.list <- validateGeneOrder(expr.x.df=exprs(expr.list[['CCLE']]),
                                  expr.y.df=exprs(expr.list[['GDSC']]))
genedist.list <- list()
genedist.list[['match']] <- generateGeneDist(expr.xy.list, type='match')
genedist.list[['mismatch']] <- generateGeneDist(expr.xy.list, type='mismatch')

#Rank the genes according to variance
gene.universe <- unlist(genedist.list[['match']])
var.universe <- (sd(gene.universe, na.rm=TRUE) ^ 2)

var.vals <- sapply(genedist.list[['match']], function(x) sd(x, na.rm=TRUE)^2)
# subset.genes <- names(var.vals)[order(var.vals)[c(1:3000)]]   # Obsolete, doesn't do what i thought at the time it does
# subset.genes.tail <- names(var.vals)[order(var.vals)[c(7000:10000)]]  # Obsolete, doesn't do what i thought at the time it does
subset.genes <- names(sort(var.vals, decreasing = TRUE))[1:3000]   # New/updated
subset.genes.tail <- names(sort(var.vals, decreasing = FALSE))[1:3000]  # New/updated


# Does an f-test between each genes variance to the gene universe variance
# 10845 genes with a qval < 0.01 out of 11351, something went wrong
{
  pval.gene.uni <- sapply(c(1:length(genedist.list[['match']])), function(x){
    ftest.pval <- var.test(genedist.list[['match']][[x]], 
                           gene.universe, alternative = 'less')$p.val
    return(ftest.pval)
  })
  pval.adj.gene.uni <- p.adjust(pval.gene.uni, method = "bonferroni")
  length(which(pval.adj.gene.uni < 0.01))
}

# Does an all genes by all genes pairwise f-test (alternative=two.sided... should be less)
# Generates a matrix of significantly different variances in a pairwise manner
{
  pval.mat <- sapply(c(1:length(genedist.list[['match']])), function(x){
    ftest.pval <- c()
    xy.ftest.pval <- sapply(c(1:length(genedist.list[['match']])), function(y){
      ftest.pval <- c(ftest.pval,
                      var.test(genedist.list[['match']][[x]], 
                               genedist.list[['match']][[y]])$p.val)
    })
    return(xy.ftest.pval)
  })
  pval.adj.mat <- matrix((p.adjust(pval.mat, method="bonferroni")), ncol=ncol(pval.mat))
  colnames(pval.adj.mat) <- names(genedist.list[['match']])
  rownames(pval.adj.mat) <- names(genedist.list[['match']])
}


# Plots the density distributions between two genes given their index
{
  plotX <- function(x,y){
    plot(density(genedist.list[['match']][[x]], na.rm=TRUE))
    lines(density(genedist.list[['match']][[y]], na.rm=TRUE), col="red")
    text(x=c(0.5,0.5), y=c(0.2, 0.1), col=c("black", "red"),
         labels=c(round(sd(genedist.list[['match']][[x]], na.rm=TRUE), 4),
                  round(sd(genedist.list[['match']][[y]], na.rm=TRUE), 4)))
  }
  
  plotX(x=7127,
        y=11194)
}
expr.list.bkup <- expr.list
drug.auc.list.bkup <- drug.auc.list
expr.list <- list()
expr.list[['CGP']] <- summarizeMolecularProfiles(GDSC, mDataType="rna",
                                                 cell.lines=cellNames(GDSC))
expr.list[['GDSC']] <- summarizeMolecularProfiles(GDSC, mDataType="rna2",
                                                  cell.lines=cellNames(GDSC))
expr.list[['CCLE']] <- summarizeMolecularProfiles(CCLE, mDataType="rna",
                                                  cell.lines=cellNames(CCLE))

# Get sensitivity profiles
drug.auc.list <- list()
drug.auc.list[['GDSC']] <- summarizeSensitivityProfiles(
  pSet=GDSC,
  sensitivity.measure='auc_recomputed',
  summary.stat='median')
drug.auc.list[['CCLE']] <- summarizeSensitivityProfiles(
  pSet=CCLE,
  sensitivity.measure='auc_recomputed',
  summary.stat='median')

ds.col <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99")
names(ds.col) <- c("GDSC", "CGP", "CCLE", "PFIZER")
outdir <- '/Users/rquevedo/Desktop/bhk_lab/results/phenotypes/drug_expr-plots'
comparePhenotypes(cl1="HPAC", ds1="CCLE", debug.mode = TRUE, loading.factors=NA, 
                  cl2="KCI-MOH1", ds2="CCLE", expr.ids=l1000.df$Entrez.Gene.ID, 
                  input.type='ENTREZID', gen.plots=FALSE, outdir=outdir,
                  rm.na=FALSE, expr.idx=subset.genes, 
                  mark.genes=seg.genes.list[['hscra1']]$Ensembl)
comparePhenotypes(cl1="HPAC", ds1="CCLE", debug.mode = TRUE, loading.factors=NA, 
                  cl2="KCI-MOH1", ds2="CCLE", expr.ids=l1000.df$Entrez.Gene.ID, 
                  input.type='ENTREZID', gen.plots=FALSE, outdir=outdir,
                  rm.na=FALSE, expr.idx=subset.genes.tail, 
                  mark.genes=seg.genes.list[['hscra2']]$Ensembl)

cl.id <- "REH" # 'Hs 940.T' # 'NCI-H23' #"REH" #"PC-3"  #"SR"
comparePhenotypes(cl1=cl.id, ds1="GDSC", debug.mode = TRUE, loading.factors=NA, 
                  cl2=cl.id, ds2="CCLE", expr.ids=l1000.df$Entrez.Gene.ID, 
                  input.type='ENTREZID', gen.plots=FALSE, outdir=outdir,
                  rm.na=FALSE, expr.idx=subset.genes)
comparePhenotypes(cl1=cl.id, ds1="GDSC", debug.mode = TRUE, loading.factors=NA, 
                  cl2=cl.id, ds2="CCLE", expr.ids=l1000.df$Entrez.Gene.ID, 
                  input.type='ENTREZID', gen.plots=FALSE, outdir=outdir,
                  rm.na=FALSE, expr.idx=subset.genes.tail, 
                  mark.genes=seg.genes.list[['hscra2']]$Ensembl)


comparePhenotypes(cl1="A253", ds1="GDSC", debug.mode = TRUE, loading.factors=NA, 
                  cl2="A253", ds2="CCLE", expr.ids=l1000.df$Entrez.Gene.ID, 
                  input.type='ENTREZID', gen.plots=FALSE, outdir=outdir,
                  rm.na=FALSE, expr.idx=subset.genes, mark.genes=seg.genes$Ensembl)
comparePhenotypes(cl1="A253", ds1="GDSC", debug.mode = TRUE, loading.factors=NA, 
                  cl2="A253", ds2="CCLE", expr.ids=l1000.df$Entrez.Gene.ID, 
                  input.type='ENTREZID', gen.plots=FALSE, outdir=outdir,
                  rm.na=FALSE, expr.idx=subset.genes.tail, mark.genes=seg.genes$Ensembl)

hscr.data <- list()
base.dir <- '/Users/rquevedo/Desktop/bhk_lab/results/phenotypes/summary_metrics/'
tsv.file <-  "GDSC-REH.CCLE-REH" #"GDSC-NCI-H23.CCLE-NCI-H23" 
#"GDSC-REH.CCLE-REH"  #GDSC-PC-3.CCLE-PC-3
hscr.data[['hscra1']] <- read.csv(file.path(base.dir, 
                                            paste("hscra1", tsv.file, "tsv", sep=".")), 
                                  sep="\t", header = TRUE)
hscr.data[['hscra2']] <- read.csv(file.path(base.dir, 
                                            paste("hscra2", tsv.file, "tsv", sep=".")), 
                                  sep="\t", header = TRUE)
exprArr.anno <- read.csv("/Users/rquevedo/git/ccl_phenotype/reference/HT_HG-U133A.na36.annot.csv/HT_HG-U133A.na36.annot.txt",
                         sep="\t", header = TRUE)
loci <- gsub(" \\([+-].+", "", exprArr.anno$Alignments, perl=TRUE)
exprArr.anno$chr <- gsub("chr", "", sapply(strsplit(loci, split=":"), function(x) x[[1]]))
exprArr.anno$start <- as.integer(gsub("chr.+:", "", sapply(strsplit(loci, split="-"), function(x) x[[1]]), perl=TRUE))
exprArr.anno$end <- as.integer(sapply(strsplit(loci, split="-"), function(x) x[[2]]))
exprArr.anno <- exprArr.anno[order(exprArr.anno$start),]
exprArr.anno.list <- split(exprArr.anno, exprArr.anno$chr)

quant.var <- 0.95
seg.genes.list <- lapply(hscr.data, function(each.hscr.df) {
  quant.thresh <- quantile(each.hscr.df$adjusted_diff, quant.var)
  
  each.hscr.list <- split(each.hscr.df, each.hscr.df$Chr)
  #1
  chr.hscr.list <- lapply(each.hscr.list, function(chr.hscr.df){
    sub.hscr.df <- chr.hscr.df[with(chr.hscr.df, which(adjusted_diff > quant.thresh)),]
    if(nrow(sub.hscr.df) > 0){
      seg.start.idx <- c(1, 1+which(diff(sub.hscr.df$adjusted_diff)!=0))
      seg.end.idx <-   c((seg.start.idx - 1)[-1], nrow(sub.hscr.df))
      
      sub.reformat.hscr.df <- data.frame("Start.bin" = sub.hscr.df[seg.start.idx,]$Start.bin,
                                         "End.bin" = sub.hscr.df[seg.end.idx,]$End.bin,
                                         "Chr" = sub.hscr.df[seg.start.idx,]$Chr,
                                         "adjusted_diff" = sub.hscr.df[seg.start.idx,]$adjusted_diff,
                                         "cn1"=sub.hscr.df[seg.start.idx,]$cn1,
                                         "cn2"=sub.hscr.df[seg.start.idx,]$cn2,
                                         "wf"=sub.hscr.df[seg.start.idx,]$wf)
    } else {
      sub.reformat.hscr.df <- sub.hscr.df
    }
    return(sub.reformat.hscr.df)
  })
  
  seg.genes.list <- lapply(names(chr.hscr.list), function(chr.id){
    seg.df <- chr.hscr.list[[as.character(chr.id)]]
    ref.df <- exprArr.anno.list[[as.character(chr.id)]]
    
    if(nrow(seg.df) > 0){
      seg.intervals <- Intervals(as.matrix(seg.df[,c('Start.bin', 'End.bin')]))
      ref.intervals <- Intervals(as.matrix(ref.df[,c("start", "end")]))
      interval_overlap(seg.intervals, ref.intervals)
      chr.seg.genes <- ref.df[unlist(interval_overlap(seg.intervals, ref.intervals)), c("Alignments", "Probe.Set.ID", "Gene.Symbol", "Ensembl")]
      chr.seg.genes$Ensembl <- gsub(" \\/.*", "", chr.seg.genes$Ensembl)
    } else {
      chr.seg.genes <- data.frame("Alignments"=NA, "Probe.Set.ID"=NA, "Gene.Symbol"=NA, "Ensembl"=NA)
    }
    return(chr.seg.genes)
  })
  
  seg.genes <- do.call("rbind", seg.genes.list)
  return(seg.genes)
  #mark.genes <- paste(seg.genes$Ensembl, "at", sep="_")
  
  #chr.hscr.df <- do.call("rbind", chr.hscr.list)
})


# # Subset the Expr data for the new set of Genes
# expr.custom.list <- list()
# expr.custom.list[['CGP']] <- exprs(expr.list[['CGP']])[match(subset.genes, rownames(exprs(expr.list[['CGP']]))),]
# expr.custom.list[['GDSC']] <- exprs(expr.list[['GDSC']])[match(subset.genes, rownames(exprs(expr.list[['GDSC']]))),]
# expr.custom.list[['CCLE']] <- exprs(expr.list[['CCLE']])[match(subset.genes, rownames(exprs(expr.list[['CCLE']]))),]


# --------------------------------------------------------------------------------

getResidual <- function(auc.x, auc.y, cell.ids, analysis.type='concordant', filt=FALSE){
  if(nrow(auc.x) == nrow(auc.y)){
    if(analysis.type == 'concordant'){
      delta.auc.list <- list()
      auc.x <- auc.x[,cell.ids]
      auc.y <- auc.y[,cell.ids]
      
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

ref.ids <- colnames(drug.auc.list[['CCLE']])[-which(colnames(drug.auc.list[['CCLE']]) %in% all.cell.ids)]

analysis.type = 'concordant'

ds1 <- 'GDSC'
ds2 <- 'CCLE'
cl1 <- ref.ids[1]
cell.ids <- ref.ids
cell.ids <- NA

subset.genes
subset.genes.tail

cl.id <- 'ChaGo-K-1'
cl.id <- 'NCI-H2126'
cl.id <- 'REH'
cl.id <- 'LoVo'
expr.mat <- getExprMat(cl1, ds1, cl2, ds2)
expr.idx <- paste(mark.genes[-which(mark.genes %in% "---" | is.na(mark.genes))], "_at", sep="")
row.idx <- match(paste(mark.genes, "at", sep="_"), rownames(expr.mat))
row.idx <- row.idx[!is.na(row.idx)]

getSOR <- function(expr.mat, expr.idx){
  expr.mat.tmp <- subsetExprMat(expr.mat=expr.mat, expr.ids=NA, 
                                input.type=NA, expr.idx=expr.idx)
  norm.expr.mat <- normalizeMedianAbsValues(expr.mat.tmp)
  dir.residuals <- getSqResiduals(norm.expr.mat, analysis='diagonal')
  return(sum(dir.residuals$residual))
}
plRes()

# Concordant:  Get the Expression matrix for cl1 and cl2 for ds1 and ds2
#  as well as obtaining the square of residuals for each of the expression matrix
if(analysis.type == 'concordant'){
  sor.list <- lapply(cell.ids, function(cl.id){
    cat(paste("Concordance:", cl.id, "\n"))
    expr.mat <- getExprMat(cl.id, ds1, cl.id, ds2)
    
    sor.ds <- tryCatch(data.frame("tail" = getSOR(expr.mat, expr.idx=subset.genes.tail),
                                  "head" = getSOR(expr.mat, expr.idx=subset.genes)), 
                       error=function(e){ NA })
    if(!is.na(sor.ds)) return(sor.ds)
  })
  
  names(sor.list) <- cell.ids
  sor.df <- do.call("rbind", sor.list)
} else if (analysis.type == 'discordant' && is.na(cell.ids)){
  delta.auc <- sapply(colnames(auc.x), function(each.cl.x){
    combine.expr.row <- sapply(colnames(auc.y), function(each.cl.y){
      if(each.cl.x != each.cl.y){
      }
    })
  })
}


# Plot Residuals Function:  Creates a scatterplot that plots the residuals according to a 
# theoretical 1:1 GDSC to CCLE expression
plRes <- function(){
  # split.screen(c(2,1))
  # screen(1)
  # plot(density(dir.residuals$residual))
  with(dir.residuals, plot(raw.x, raw.y, xlim=c(0,14), ylim=c(0,14), main=cl.id))
  with(dir.residuals, points(xval, yval, col="red", pch=16))
  with(dir.residuals, segments(x0=xval, y0=yval, 
                               x1=raw.x, y1=raw.y, col=alpha("red", 0.5)))
  # close.screen(all.screens=TRUE)
}

# mark.genes=seg.genes.list[['hscra1']]$Ensembl
# row.idx <- match(paste(mark.genes, "at", sep="_"), rownames(expr.mat))
# row.idx <- row.idx[!is.na(row.idx)]
# expr.mat <- expr.mat[row.idx,]
dir.residuals <- getSqResiduals(expr.top.mat, analysis='diagonal')
sum(dir.residuals$residual)
dir.residuals <- getSqResiduals(expr.bot.mat, analysis='diagonal')
sum.of.residuals <- sum(dir.residuals$residual)
sor.vals <- c(sor.vals, sum.of.residuals)






comparePhenotypes(cl1="HPAC", ds1="CCLE", debug.mode = TRUE, loading.factors=NA, 
                  cl2="KCI-MOH1", ds2="CCLE", expr.ids=l1000.df$Entrez.Gene.ID, 
                  input.type='ENTREZID', gen.plots=FALSE, outdir=outdir,
                  rm.na=FALSE, expr.idx=subset.genes, 
                  mark.genes=seg.genes.list[['hscra1']]$Ensembl)


cnv.ids <- unique(sort(unlist(lapply(cnvDiff.list, names))))
cnv.ids <- c("HT-29", "HuP-T4", "MCF7", "Mewo", "NCI-H23", "REH", "UACC-812")
cnvDiff.delta.auc <- getDeltaAuc(auc.x=drug.auc.list[['CCLE']], 
                                 auc.y=drug.auc.list[['GDSC']], 
                                 cell.ids=cnv.ids, analysis.type='concordant',
                                 filt=FALSE)


# --------------------------------------------------------------------------------

load("/Users/rquevedo/Desktop/bhk_lab/data/cnv_probeset/copy_ratio.probeset.RData")

cr.mat.summ <- round(cr.mat.summ, 4)

# Probeset priority: http://www.affymetrix.com/support/help/faqs/mouse_430/faq_8.jsp
# priority for duplicates _at, _s_at, _x_at
dup.ensg.id <- which(table(ord.df.summ$ensg.id) > 1, useNames = TRUE) # Identify ENSG ids that are duplicated
dup.ensg.df <- lapply(names(dup.ensg.id), function(x){
  pattern.ord <- list("a" = "[0-9]_at",
                      "s" = "[0-9]_s_at",
                      "x" = "[0-9]_x_at")
  dup.ensg.df <- ord.df.summ[which(ord.df.summ$ensg.id == x),]
  asx.order <- lapply(pattern.ord, function(pat) grep(pat, dup.ensg.df$Probe.Set.ID)[1])
  
  new.ensg.df <- dup.ensg.df[asx.order[['a']],]
  if(is.na(asx.order[['a']])) new.ensg.df <- dup.ensg.df[asx.order[['s']],]
  if(is.na(asx.order[['a']]) && is.na(asx.order[['s']])) new.ensg.df <- dup.ensg.df[asx.order[['x']],]
  if(all(sapply(asx.order, is.na)))  new.ensg.df <- dup.ensg.df[1,] #print(paste("Warning: No match for", x))
  return(new.ensg.df)
}) # Of the duplicate ENSG, prioritize based on probeset specificity
dup.ensg.df <- do.call("rbind", dup.ensg.df)

# Obtain the "row.ids" for all initial and duplicated ENSG ID values
duplicated.rowids <- sapply(names(dup.ensg.id), function(x){
  return(rownames(ord.df.summ[which(ord.df.summ$ensg.id == x),]))
})

orig.order.ord.df <- rownames(ord.df.summ)
ord.df.summ <- ord.df.summ[-match(unlist(duplicated.rowids), rownames(ord.df.summ)),] # Remove all duplicate ENSGid rows
ord.df.summ <- rbind(ord.df.summ, dup.ensg.df) # Append the prioritized probeset ENSG ids
ord.df.summ <- ord.df.summ[-which(is.na(ord.df.summ$ensg.id)),] # Remove NA ENSG Ids


# Intersect all ENSG IDs between the two expression sets and the annotated CNV sets
ds1.ensg.ids <- rownames(exprs(expr.list[[ds1]]))
ds2.ensg.ids <- rownames(exprs(expr.list[[ds2]]))
all.ensgid.interesct <- sort(Reduce(intersect, list(ds1.ensg.ids, 
                                                    ds2.ensg.ids, 
                                                    ord.df.summ$ensg.id)))

# Subset based on the ENSG Ids shared between all exprs and cnv sets
ds1.exprs.df <- exprs(expr.list[[ds1]])[all.ensgid.interesct,]
ds2.exprs.df <- exprs(expr.list[[ds2]])[all.ensgid.interesct,]
ord.df.summ <- ord.df.summ[match(all.ensgid.interesct, ord.df.summ$ensg.id),]
cr.mat.summ <- cr.mat.summ[match(rownames(ord.df.summ), orig.order.ord.df),]

# Load cell.line.anno and subset for just a specific tissue
load("/Users/rquevedo/git/cnv_fingerprint/v2.0/data/merged_annotations.Rdata")
tissue.cell.line.anno <- cell.line.anno[which(cell.line.anno$unique.tissueid %in% 'haematopoietic_and_lymphoid_tissue'),]
getCclTissueIds <- function(x, tissue.cell.line.anno, col.id='unique.cellid'){
  if(col.id != 'unique.cellid'){
    # For the cnv sets which uses colnames of 'CEL_filename.hscr.a1.RData'
    colnames(x) <- gsub(".hscr.*$", ".CEL", colnames(x), ignore.case = TRUE)
    expr.celids <- sapply(colnames(x), function(y) which(y == tissue.cell.line.anno, 
                                                         arr.ind = TRUE))  # Match .CEL Files
    id.df <- as.data.frame(do.call("rbind", expr.celids), row.names = FALSE)
    id.df$unique.cellid <- tissue.cell.line.anno[id.df$row, 'unique.cellid']
    id.df$dataset <- colnames(tissue.cell.line.anno)[id.df$col]
    id.df$cel.id <- sapply(c(1:nrow(id.df)), function(y){
      gsub(".CEL", "", tissue.cell.line.anno[id.df[y, 'row'], id.df[y, 'dataset']], 
           ignore.case=TRUE)
    })
    id.df <- id.df[which(id.df$col %in% col.id),] #Subset for just DS1 or DS2 filenames (e.g. excludes pfizer)
    
    # id.df dataframe contains the index in cell.line.anno, as well as dataset, unique.cellid, and cel filename
    uniq.ids <- id.df[,c('cel.id', 'unique.cellid')]
  } else {
    # For the exprs sets which uses colnames of 'unique.cellid'
    expr.celids <- match(colnames(x), unlist(tissue.cell.line.anno[,col.id]))
    uniq.ids <- tissue.cell.line.anno[expr.celids[!is.na(expr.celids)], col.id]
  } 
  
  return(uniq.ids)
}

# Obtain a vector of all CCL identifiers in datasets that match the tissue of interest
ds1.expr.celids <- getCclTissueIds(ds1.exprs.df, tissue.cell.line.anno)
ds2.expr.celids <- getCclTissueIds(ds2.exprs.df, tissue.cell.line.anno)
cn.expr.celids <- getCclTissueIds(cr.mat.summ, 
                                  tissue.cell.line.anno[match(intersect(ds1.expr.celids, ds2.expr.celids), 
                                                              tissue.cell.line.anno$unique.cellid),], 
                                  sapply(c(ds1, ds2), function(y) {
                                    grep(paste(y, "filename", sep="."), 
                                         colnames(tissue.cell.line.anno), ignore.case=TRUE)
                                  }))
# Find CCLs that are found between all exprs and cnv datasets
expr.cn.ccl.overlap <- Reduce(intersect, list(ds1.expr.celids, ds2.expr.celids, cn.expr.celids$unique.cellid))


# Plots copy.ratio values against expression for all probesets in a single cell line
pdf("/Users/rquevedo/Desktop/expr.cn.cclPlots.highCorr.pdf")
lapply(expr.cn.ccl.overlap, function(x){
  cel.id <- cn.expr.celids[match(x, cn.expr.celids$unique.cellid), 'cel.id']
  col.idx <- grep(paste("^", cel.id, ".hscr.*", sep=""), colnames(cr.mat.summ), ignore.case=TRUE)
  cr.expr.df <- data.frame("copy.ratio"=cr.mat.summ[,col.idx], 
                           "ds1.expr"=ds1.exprs.df[,x])
  if(length(col.idx) > 0){
    m <- ggplot(cr.expr.df, aes(x = copy.ratio, y = ds1.expr)) +
      xlim(0, 2) + ylim(0, 15) + ggtitle(x)
    m + geom_point() + stat_density2d(aes(fill = ..level..), geom="polygon")
  }
})
dev.off()


# Plots copy.ratio values against expression for one probeset across all tissue-specific cell lines
pdf("/Users/rquevedo/Desktop/bhk_lab/results/phenotypes/summary_metrics/expr.cn.probesetPlots.pdf")
min.cor <- 0.75
probesets.per.plot <- 50
y.val <- 15
plot.cnt <- 1
cn.sensitive.expr.df <- data.frame()
for(each.probeset in c(1:nrow(ds1.exprs.df))){
  probeset.expr.cn.df <- lapply(expr.cn.ccl.overlap, function(x){
    cel.id <- cn.expr.celids[match(x, cn.expr.celids$unique.cellid), 'cel.id']
    col.idx <- grep(paste("^", cel.id, ".hscr.*", sep=""), colnames(cr.mat.summ), ignore.case=TRUE)
    if(length(col.idx) == 0) col.idx <- 1
    return(data.frame("copy.ratio"=cr.mat.summ[each.probeset, col.idx], 
                      "ds1.expr"=ds1.exprs.df[each.probeset,x]))
  })
  probeset.expr.cn.df <- do.call("rbind", probeset.expr.cn.df)
  
  
  max.D <- round(max(cooks.distance(glm(probeset.expr.cn.df))),2)
  R.val <- round(cor(probeset.expr.cn.df, use="complete.obs")[1,2], 2)
  cn.sensitive.expr.df <- rbind(cn.sensitive.expr.df,
                                data.frame("probeset"=rownames(ds1.exprs.df)[each.probeset],
                                           "r.val"=R.val,
                                           "max.d"=max.D))
  if(cor(probeset.expr.cn.df, use="complete.obs")[1,2] > min.cor){
    probeset.col <- rainbow(100)[plot.cnt]
    text(x = 0, y=y.val, adj=1, pos = 4, col=probeset.col, cex = 0.6,
         labels = paste(rownames(ds1.exprs.df)[each.probeset], 
                        ", R: ", R.val,
                        ", max.D: ", max.D, sep=""))
    y.val <- y.val - 0.5
  } else {
    probeset.col <- alpha("grey", 0.2)
    #probeset.expr.cn.df <- data.frame()
  }
  
  if(plot.cnt == 1) {
    plot(probeset.expr.cn.df, col=probeset.col,
         xlim=c(0,2), ylim=c(0,15))
    with(probeset.expr.cn.df,
         abline(lm(ds1.expr ~ copy.ratio), col=probeset.col))
    plot.cnt <- plot.cnt + 1
  } else if (plot.cnt > (probesets.per.plot - 1)) {
    points(probeset.expr.cn.df, col=probeset.col)
    with(probeset.expr.cn.df,
         abline(lm(ds1.expr ~ copy.ratio), col=probeset.col))
    plot.cnt <- 1
    y.val <- 15
  } else {
    points(probeset.expr.cn.df, col=probeset.col)
    with(probeset.expr.cn.df,
         abline(lm(ds1.expr ~ copy.ratio), col=probeset.col))
    plot.cnt <- plot.cnt + 1
  }
}
dev.off()

# Get probeset IDs that are above a 95th quantile of correlation
cn.sensitive.expr.df <- cn.sensitive.expr.df[order(cn.sensitive.expr.df$r.val, decreasing = FALSE),]
r.val.idx <- cn.sensitive.expr.df$r.val > quantile(cn.sensitive.expr.df$r.val, 0.90)
highR.cnv.expr.df <- cn.sensitive.expr.df[which(r.val.idx), 'probeset']  

# Plots copy.ratio values and expression against genomic loci, intensity colour scaled to r.val
load("/Users/rquevedo/git/reference/cytoband/chr_cytoband.hg19.Rdata")  #loading chr.df to plot genomic loci
pdf("/Users/rquevedo/Desktop/expr.cn.highRprobeset.lociPlots.pdf") #highRprobeset, rankedR
for(each.ccl in expr.cn.ccl.overlap){
  expr.col.idx <- cn.expr.celids[match(each.ccl, cn.expr.celids$unique.cellid), 'cel.id']
  cnv.col.idx <- grep(paste("^", expr.col.idx, ".hscr.*", sep=""), colnames(cr.mat.summ), ignore.case=TRUE)
  cnv.col <- 'red'
  expr.col <- 'blue'
  
  plot(0, type='n', ylim=c(0,15), xlim=c(0, chr.df$cumEnd[22]), main=each.ccl)
  abline(v=chr.df$cumStart, col="grey")
  
  cbs.df <- data.frame()
  #for(each.probeset in cn.sensitive.expr.df$probeset){
  for(each.probeset in highR.cnv.expr.df){
    ensg.idx <- match(each.probeset, ord.df.summ$ensg.id)
    r.val <- abs(cn.sensitive.expr.df[which(cn.sensitive.expr.df$probeset == each.probeset),]$r.val)
    
    chr.df.row <- chr.df[match(ord.df.summ[ensg.idx, 'Chr'], chr.df$chr),]
    probe.start <- ord.df.summ[ensg.idx,]$Start.bin + chr.df.row$cumStart
    probe.end <- ord.df.summ[ensg.idx,]$End.bin + chr.df.row$cumStart
    
    cr.val <- cr.mat.summ[ensg.idx, cnv.col.idx]
    expr.val <- ds1.exprs.df[each.probeset, each.ccl]
    
    segments(x0=probe.start, y0=cr.val, 
             x1=probe.end, y1=cr.val, col=alpha(cnv.col, r.val))
    segments(x0=probe.start, y0=expr.val, 
             x1=probe.end, y1=expr.val, col=alpha(expr.col, r.val))
    
    cbs.df <- rbind(cbs.df, data.frame("probe.loc"=probe.start,
                                       "cr.val"=cr.val,
                                       "expr.val"=expr.val))
  }
  
  for(i in c("cr", "expr")){
    if (i == 'cr') {
      cbs.col <- cnv.col
      in.val <- "cr.val"
      a.val <- 0.05
    } else {
      cbs.col <- expr.col
      in.val <- "expr.val"
      a.val <- 0.25
    }

    cbs.out.df <- tryCatch({
      cbs.out <- getCbsSeg(chr=as.character(1), val=cbs.df[,in.val],
                           pos = cbs.df$probe.loc, name="test", aval=a.val)
      cbs.out.df <- cbs.out$output
      
      for(i in seq(1:dim(cbs.out.df)[1])) {
        #guideline at CNV median
        segments(cbs.out.df[i,'loc.start'],
                 as.numeric(cbs.out.df[i,'seg.mean']),
                 cbs.out.df[i,'loc.end'],
                 as.numeric(cbs.out.df[i,'seg.mean']),
                 col=cbs.col,
                 lwd=2)       
      }          
    }, error=function(e){
      return(NA)
    })
  }
}
dev.off()




# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# cc <- cellNames(common[[1]])
# ge.cor <- sapply(cc, function (x, d1, d2) {
#   return (cor(d1[ ,x], d2[ ,x], method="pearson",
#                      use="pairwise.complete.obs"))
# }, d1=expr.sub.list[['CGP']], d2=expr.sub.list[['CCLE']])
# ge.mismatch.cor <- sapply(cc, function(x) {
#   sapply(cc, function(y, d1, d2) {
#     if(!x %in% y){
#       return (stats::cor(d1[ , x], d2[ , y], method="pearson",
#                          use="pairwise.complete.obs"))
#     } else {
#       return(1)
#     }
#   }, d1=expr.sub.list[['CGP']], d2=expr.sub.list[['CCLE']])
# })
# 
# box.info <- boxplot(data.frame("match"=ge.cor,
#                                "mismatch"=as.vector(ge.mismatch.cor)), 
#                     ylim=c(0,1), main="GDSC-CCLE Correlation", 
#                     sub="GDSC:rna; CCLE:rna", ylab="Spearman Correlation")
