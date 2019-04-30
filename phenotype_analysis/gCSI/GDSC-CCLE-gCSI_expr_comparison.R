packages <- c("AnnotationDbi", "org.Hs.eg.db",
              "PharmacoGx", "scales")
tmp <- suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE))

outdir <- '~/Desktop/bhk_lab/results/phenotypes/expr_plots'
refdir <- '~/git/reference/l1000'

outdir <- '/results/pharmaco'  #codeocean
refdir <- '/data/ref'  #codeocean

dir.create(outdir, recursive=TRUE)
load(file.path(refdir, "l1000.Rdata"))

#### FUNCTIONS ####
plotHist <- function(r, Y, lower.y=0.001, upper.y=0.005){
  h <- hist(r, breaks=50, plot=FALSE)
  cuts <- cut(h$breaks, c(-Inf, (Y-lower.y), (Y + upper.y), Inf))
  plot(h , col=c("grey", "red", "grey")[cuts], xlim=c(0,1),
       xlab="Pearson correlation", main="GDSC-CCLE L1000 genes")
  text(x=0.2, y=quantile(h$counts, 0.95),
       labels=paste0("r = ",round(Y, 3)))
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



## Takes a correlation between all matching cell lines (grep) between all pairwise datasets
corrExprMatching <- function(expr, ds.ids){
    sapply(ds.ids, function(a){
        sapply(ds.ids, function(b){
            sapply(colnames(expr[[a]]), function(clA){
                clB <- grep(paste0("^", clA, "$"), colnames(expr[[b]]))
                cor(expr[[a]][,clA], expr[[b]][,clB],
                    use = "na.or.complete", method = "pearson")
            })
        })
    })
}

## Summarizes (mean and sd) of the pairwise correlations
summarizeCorr <- function(corr){
    mean.mat <- sapply(corr, function(dsA){
        apply(dsA, 2, function(dsB){
            dsB <- suppressWarnings(as.numeric(as.character(dsB)))
            "mean"=round(mean(dsB, na.rm=TRUE), 3)
        })
    })

    sd.mat <- sapply(corr, function(dsA){
        apply(dsA, 2, function(dsB){
            dsB <- suppressWarnings(as.numeric(as.character(dsB)))
            "sd"=round(sd(dsB, na.rm=TRUE),3)
        })
    })
    return(list("mean"=mean.mat,
               "sd"=sd.mat))
}

#### MAIN ####
#### Obtain PharmacoGx Data ####
#Download PSets
#availablePSets()
GDSC1000 <- downloadPSet("GDSC1000")
gCSI <- downloadPSet("gCSI")
CCLE <-downloadPSet("CCLE")

#   PharmacoSets for cell lines found in both GDSC, CCLE, and gCSI
commonCl <- list()
commonCl[['CCLE-GDSC']] <-  intersectPSet(list('CCLE'=CCLE,
                                               'GDSC'=GDSC1000),
                                          intersectOn=c("cell.lines"))
commonCl[['CCLE-gCSI']] <-  intersectPSet(list('CCLE'=CCLE,
                                               'gCSI'=gCSI),
                                          intersectOn=c("cell.lines"))
commonCl[['GDSC-gCSI']] <-  intersectPSet(list('GDSC'=GDSC1000,
                                               'gCSI'=gCSI),
                                          intersectOn=c("cell.lines"))
commonCl[['all']] <- intersectPSet(list('CCLE'=CCLE,
                                        'GDSC'=GDSC1000,
                                        "gCSI"=gCSI),
                                   intersectOn=c("cell.lines"))

# Find all common genes between rna datasets
commonGenes <- list()
commonGenes[['CCLE-GDSC']] <- intersect(fNames(CCLE, "rna"),
                                        fNames(GDSC1000,"rna"))
commonGenes[['CCLE-gCSI']] <- intersect(fNames(CCLE, "rna"),
                                        fNames(gCSI,"rnaseq"))
commonGenes[['GDSC-gCSI']] <- intersect(fNames(GDSC1000, "rna"),
                                        fNames(gCSI,"rnaseq"))
commonGenes[['all']] <- Reduce(function(d1, d2) intersect(d1, d2),
    list(fNames(GDSC1000, "rna"),
         fNames(CCLE,"rna"),
         fNames(gCSI,"rnaseq")))

print("Number of common genes in RNAseq, per dataset: ")
print(sapply(commonGenes, length))


# Subsetting the RNA datasets just for the genes found in all datasets
expr.list <- list()
expr.list[['CCLE']] <- getExprVals(CCLE, NULL, "rna", commonGenes[['all']])
expr.list[['GDSC']] <- getExprVals(GDSC1000, NULL, "rna", commonGenes[['all']])
expr.list[['gCSI']] <- getExprVals(gCSI, NULL, "rnaseq", commonGenes[['all']])


# Find the L1000 Genes Expression values for overlapping cell lines:
expr.sub.list <- list()
expr.sub.list[['GDSC']] <- subsetExprMat(exprs(expr.list[['GDSC']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)
expr.sub.list[['CCLE']] <- subsetExprMat(exprs(expr.list[['CCLE']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)
expr.sub.list[['gCSI']] <- subsetExprMat(exprs(expr.list[['gCSI']]), l1000.df$Entrez.Gene.ID, input.type='ENTREZID', rm.na=FALSE)


expr.list <- lapply(expr.list, function(x) exprs(x))


#### Generating correlations ####
# Pearson correlation of L1000 genes
ds.ids <- names(expr.sub.list)

l1000.corr <- corrExprMatching(expr.sub.list, ds.ids)
all.corr <- corrExprMatching(expr.list, ds.ids)

# Correlations between gene expression of L1000 and All genes
print("Correlation values of L1000 gene expression: ")
print(summarizeCorr(l1000.corr))
print("Correlation values of All overlapping genes expression: ")
print(summarizeCorr(all.corr))



#### Visualization: Correlation histogram ####
# Pearson correlation of L1000 genes for tissue-specificy
l1000.Y <- cor(expr.sub.list[['GDSC']][,'HPAC'],
               expr.sub.list[['CCLE']][,'KCI-MOH1'],
               use = "na.or.complete", method = "pearson")

all.Y <- cor(expr.list[['GDSC']][,'HPAC'],
             expr.list[['CCLE']][,'KCI-MOH1'],
             use = "na.or.complete", method = "pearson")

pdf(file.path(outdir, "kci-moh1.l1000.pdf"))
plotHist(unlist(l1000.corr[['GDSC']][,'CCLE']), l1000.Y)
dev.off()

pdf(file.path(outdir, "kci-moh1.all.pdf"))
plotHist(unlist(all.corr[['GDSC']][,'CCLE']), all.Y)
dev.off()

# Pearson correlation of L1000 genes for tissue-specificy
l1000.Y <- cor(expr.sub.list[['GDSC']][,'KARPAS-422'],
               expr.sub.list[['CCLE']][,'OCI-LY10'],
               use = "na.or.complete", method = "pearson")

all.Y <- cor(expr.list[['GDSC']][,'KARPAS-422'],
             expr.list[['CCLE']][,'OCI-LY10'],
             use = "na.or.complete", method = "pearson")

pdf(file.path(outdir, "karpas-oci.all.pdf"))
plotHist(unlist(all.corr[['GDSC']][,'CCLE']), all.Y, upper.y=0.002)
dev.off()


