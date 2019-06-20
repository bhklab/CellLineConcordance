library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
genes.txdb = genes(txdb)

###################
#### Variables ####
PDIR='/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/seg_files'
setwd(PDIR)

load("merged_annotations.Rdata")
gcsi.metadata <- '/mnt/work1/users/bhklab/Data/GNE/genentech.snparray.manifest.csv'
anno.name <- 'gCSI'


#
seg.dir <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/seg_files'
file <- 
anno.name <- 'CCLE'
gcsi.exprsFile <- file.path(seg.dir, file)
gcsi.eset.file <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2017_GNE/PSet/gCSI.R'
source("~/git/CellLineConcordance/eSet/R/annotateGenes.R")

###################
#### Functions ####
getMapping <- function(key='ENTREZID', column=c("SYMBOL","ENSEMBL")){
  keyentrez <- keys(org.Hs.eg.db, keytype=key) 
  key.map <- select(org.Hs.eg.db, 
                    keys=keyentrez, 
                    columns = column,
                    keytype=key)
  key.map
}

# https://support.bioconductor.org/p/96427/
annotateCnvs <- function(cnv, txdb, anno=NULL, 
                         cols=c("seg.mean", "nA", "nB")){
  stopifnot(is(cnv, "GRanges"), is(txdb, "TxDb"))
  
  ## Assign EntrezID to each segment
  if(is.null(anno)) anno = genes(txdb)
  olaps = findOverlaps(cnv, anno)
  mcols(olaps)$gene_id = anno$gene_id[subjectHits(olaps)]  # Fixed the code here
  cnv_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
  cnv$gene_id = IRanges::splitAsList(mcols(olaps)$gene_id, cnv_factor)
  
  ## Assign EntrezID to each segment
  seg.entrez <- apply(as.data.frame(mcols(cnv)), 1, function(i){
    ids <- unlist(strsplit(x = as.character(unlist(i[['gene_id']])), split=","))
    segs <- do.call(rbind, replicate(length(ids), round(unlist(i[cols]),3), simplify = FALSE))
    
    as.data.frame(cbind(segs, 'ENTREZ'=ids))
  })
  seg.entrez <- do.call(rbind, seg.entrez)
  if(any(duplicated(seg.entrez$ENTREZ))) seg.entrez <- seg.entrez[-which(duplicated(seg.entrez$ENTREZ)),]
  
  
  ## Map ensembl and HUGO IDs to the ENTREZ ids
  seg.anno <- merge(seg.entrez, getMapping(), 
                    by.x="ENTREZ", by.y="ENTREZID", all.x=TRUE)
  seg.anno <- seg.anno[-which(duplicated(seg.anno$ENTREZ)),]
  for(each.col in cols){
    seg.anno[,each.col] <- as.numeric(as.character(seg.anno[,each.col])) 
  }
  
  
  list("seg"=cnv, "genes"=seg.anno)
}     

reduceEsetMats <- function(gene.lrr, cols, features='SYMBOL'){
  lapply(cols, function(each.col, features){
    print(each.col)
    keys <- c("ENTREZ", "SYMBOL", "ENSEMBL")
    m <- suppressWarnings(Reduce(f=function(x,y) merge(x,y,by=keys), 
                                 lapply(gene.lrr, function(i) i[['genes']][,c(keys, each.col)])))
    if(any(duplicated(m[,features]))) m <- m[-which(duplicated(m[,features])),]
    if(any(is.na(m[,features]))) m <- m[-which(is.na(m[,features])),]
    rownames(m) <- m[,features]
    m <- m[,-c(1:3)]
    colnames(m) <- names(cnseg.list)
    as.matrix(m)
  }, features=features)
}

##############
#### Main ####

handle <- 'CCLE'
switch(handle,
       'gCSI'={
         file <- c("gCSI.seg", "gCSI_cn.raw.tsv")
         anno.name <- 'gCSI'
         gcsi.metadata <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/seg_files/ref/genentech.snparray.manifest.csv'
         cols <- c("seg.mean", "nAraw", "nBraw")
         ids <- c('ID', 'chrom')
       },
       'CCLE'={
         file <- 'CCLE.segmaf.seg'
         anno.name <- 'CCLE'
         ccle.metadata <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/seg_files/ref/CCLE_sample_info_file_2012-10-18.txt'
         cols <- c("copy.ratio", "hscr.a1", "hscr.a2", 'modal.a1', 'modal.a2')
         ids <- c('sample', "Chromosome")
       },
       'GDSC'={
         file <- 'GDSC.segmaf.seg'
         anno.name <- 'GDSC'
         ccle.metadata <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/seg_files/ref/CCLE_sample_info_file_2012-10-18.txt'
         cols <- c("copy.ratio", "hscr.a1", "hscr.a2", 'modal.a1', 'modal.a2')
         ids <- c('sample', "Chromosome")
       })

seg <- read.table(file[1], header=TRUE, sep="\t", 
                  stringsAsFactors = FALSE, check.names=FALSE)
if(handle == 'gCSI'){
  gcsi.dat <- read.table(file[2], header=TRUE, sep="\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  gcsi.seg <- merge(seg, gcsi.dat, 
                    by.x=c("ID", "chrom", "loc.start", "loc.end"),
                    by.y=c("sample", "chr", "startpos", "endpos"))
  seg <- gcsi.seg
}

cnseg.list <- split(seg, f=seg[,ids[1]])
gene.lrr <- lapply(cnseg.list, function(each.cnseg){
  # Remove NA's and convert to a GRanges Object
  if(any(is.na(each.cnseg[,ids[2]]))) each.cnseg <- each.cnseg[-which(is.na(each.cnseg[,ids[2]])),]
  cnv <- makeGRangesFromDataFrame(each.cnseg, keep.extra.columns=TRUE, 
                                  start.field = c('loc.start', 'Start.bp'),
                                  end.field = c('loc.end', 'End.bp'))
  suppressWarnings(seqlevelsStyle(cnv) <- "UCSC")
  
  # Annotate and return the gene data.frame
  print(paste0("gCSI: ", unique(each.cnseg[,ids[1]]), " - ",
               grep(unique(each.cnseg[,ids[1]]), names(cnseg.list)), "/", length(cnseg.list)))
  anno <- suppressMessages(annotateCnvs(cnv, txdb, anno=genes.txdb, cols=cols))
  anno # list containing 'seg' (GRanges) and 'genes' (DF)
})

mats <- reduceEsetMats(gene.lrr, cols, features='SYMBOL')
save(mats, file=file.path('esets', "raw_mats", paste0(anno.name, '.Rdata')))

## Format the phenoData
if(handle == 'gCSI'){
  meta <- read.table(gcsi.metadata, sep=",", header=TRUE, 
                     stringsAsFactors = FALSE, check.names = FALSE)
  meta$INVENTORY_SAMPLE_NAME <- gsub("Jurkat,", "Jurkat", meta$INVENTORY_SAMPLE_NAME)
  meta <- as.data.frame(meta[match(names(cnseg.list), meta$INVENTORY_SAMPLE_NAME),])
  rownames(meta) <- as.character(meta$INVENTORY_SAMPLE_NAME)
} else if (handle == 'CCLE'){
  meta <- read.table(ccle.metadata, sep="\t", header=TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
  m.ids <- match(colnames(mats[[1]]), meta[,'SNP arrays'])
  meta <- meta[m.ids,]
  mats <- lapply(mats, function(i){
    non.na.idx <- which(!is.na(meta[,'SNP arrays']))
    colnames(i)[non.na.idx] <- meta[non.na.idx,'Cell line primary name']
    i
  })
  rownames(meta) <- colnames(mats[[1]])
} else if (handle == 'GDSC'){
  common.ids <- c('COSMIC.tissueid', 'unique.cellid', 'unique.tissueid')
  m.ids <- match(common.ids, colnames(cell.line.anno))
  a.ids <- grep(anno.name, colnames(cell.line.anno), ignore.case = TRUE)
  meta <- cell.line.anno[,c(m.ids, a.ids)]
  if(handle=='GDSC') meta[947,'GDSC.cellid'] <- 'TT'
  
  idx <- grep(paste0(anno.name, '.*cellid'), colnames(meta), ignore.case = TRUE)
  meta <- meta[match(colnames(mats[[1]]), meta[,idx]),]
  mats <- lapply(mats, function(i){
    non.na.idx <- which(!is.na(meta[,'unique.cellid']))
    colnames(i)[non.na.idx] <- meta[non.na.idx,'unique.cellid']
    i
  })
  rownames(meta) <- colnames(mats[[1]])
}

## Assemble the assayData environment
eset.env <- new.env()
assign("exprs", mats[[1]], envir=eset.env)
assign("A", mats[[2]], envir=eset.env)
assign("B", mats[[3]], envir=eset.env)
if(handle == 'GDSC'| handle == 'CCLE'){
  assign("modalA", mats[[4]], envir=eset.env)
  assign("modalB", mats[[5]], envir=eset.env)
}


## Assemble the eset
cl.phenoData <- new("AnnotatedDataFrame", data=meta)
cl.eset <- ExpressionSet(assayData=eset.env,
                         phenoData=cl.phenoData,
                         annotation=anno.name)
save(cl.eset, file=file.path("esets", paste0(anno.name, "_eset.Rdata")))