library("Biobase")
###################
#### Variables ####
gcsi.exprsFile <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2017_GNE/ascat'
gcsi.metadata <- '/mnt/work1/users/bhklab/Data/GNE/genentech.snparray.manifest.csv'
gcsi.eset.file <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2017_GNE/PSet/gCSI.R'
source("~/git/CellLineConcordance/eSet/R/annotateGenes.R")

###################
#### Functions ####
splitGeneLRR <- function(i){
  em <- elementMetadata(i)@listData
  seg.gene.df <- sapply(1:length(i), function(idx){
    genes <- unlist(em$gene_ids[[idx]])
    if(length(genes) > 0){
      s.gene <- cbind(rep(em$seg.mean[[idx]], length(genes)),
                      as.matrix(genes)) 
      if(any(is.na(s.gene[,2]))) s.gene <- s.gene[-which(is.na(s.gene[,2])),, drop=FALSE]
      s.gene
    }
  })
  
  seg.gene.df <- as.data.frame(do.call(rbind, seg.gene.df))
  colnames(seg.gene.df) <- c("seg.mean", "gene")
  seg.gene.df$seg.mean <- round(as.numeric(as.character(seg.gene.df$seg.mean)),3)
  seg.gene.df
}

##############
#### Main ####
summ.dir <- file.path(gcsi.exprsFile, "output", "summ")

cl.exprs <- read.table(file.path(summ.dir, "gCSI.seg"), 
                       header = TRUE, sep = "\t", 
                       stringsAsFactors = FALSE, check.names=FALSE)
cl.ids <- split(cl.exprs, f=cl.exprs$ID)
hg19.genes <- getGenes()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

## Annotate the CN segments
anno.ids <- lapply(cl.ids, function(cl.i){
  print(unique(cl.i$ID))
  cl.i <- cl.i[-which(is.na(cl.i$chrom)),]
  cl.anno <- suppressMessages(annotateSegments(cl.i, hg19.genes))
  gene.lrr <- splitGeneLRR(cl.anno)
  
  list("gene.lrr"=gene.lrr, "gr"=cl.anno)
})

## Format the exprs() matrix
cl.order.exprs <-Reduce(f = function(x,y) merge(x,y, by='gene'), 
                        lapply(anno.ids, function(i) i[['gene.lrr']]))
rownames(cl.order.exprs) <- cl.order.exprs[,1]
cl.order.exprs <- cl.order.exprs[,-1]
colnames(cl.order.exprs) <- names(anno.ids)
cl.order.exprs <- as.matrix(cl.order.exprs)

## Format the phenoData
meta <- read.table(gcsi.metadata, sep=",", header=TRUE, 
                   stringsAsFactors = FALSE, check.names = FALSE)
meta$INVENTORY_SAMPLE_NAME <- gsub("Jurkat,", "Jurkat", meta$INVENTORY_SAMPLE_NAME)
meta <- as.data.frame(meta[match(names(anno.ids), meta$INVENTORY_SAMPLE_NAME),])
rownames(meta) <- as.character(meta$INVENTORY_SAMPLE_NAME)

## Assemble the eset
anno.name <- 'gCSI'
cl.phenoData <- new("AnnotatedDataFrame", data=meta)
gcsi.eset <- ExpressionSet(assayData=cl.order.exprs,
                           phenoData=cl.phenoData,
                           annotation=anno.name)
save(gcsi.eset, file=gcsi.eset.file)
