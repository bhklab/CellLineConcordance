library("Biobase")
library(org.Hs.eg.db)
library(parallel)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

###################
#### Variables ####
handle <- 'UHN'
PDIR='/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/seg_files'
setwd(PDIR)
load('merged_annotations.Rdata')

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb.genes <- genes(txdb)

###################
#### Functions ####
getMapping <- function(in.col='ENTREZID', 
                       out.cols=c("SYMBOL", "ENSEMBL")){
  gene.map <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db, in.col), 
                     keytype="ENTREZID", columns=out.cols)
  gene.map
}

genWindowedBed <- function(bin.size=1000000, seq.style="UCSC"){
  chrs <- org.Hs.egCHRLENGTHS[c(1:22,"X", "Y")]
  
  ## Construct intervals across the genome of a certain bin size
  start.points <- seq(1, 500000000, by=bin.size)
  grl <- lapply(names(chrs), function(chr.id){
    chr <- chrs[chr.id]
    ir <- IRanges(start=start.points[start.points < chr], width=bin.size)
    end(ir[length(ir),]) <- chr
    gr <- GRanges(seqnames = chr.id, ir)
    gr
  })
  
  ## Assemble all GRanges and set seq level style
  grl <- as(grl, "GRangesList")
  suppressWarnings(seqlevelsStyle(grl) <- seq.style)
  gr <- unlist(grl)
  gr
}

segmentCNVs <- function(cnv, bed){
  olaps = findOverlaps(cnv, bed)
  
  # Flag BED bins that map to multiple CNVs
  dup.idx <- which(duplicated(subjectHits(olaps), fromLast=TRUE))
  dup.idx <- c(dup.idx, which(duplicated(subjectHits(olaps), fromLast=FALSE)))
  dup.idx <- sort(dup.idx)
  
  # Use a summary metric (Default=mean) to reduce the CNV information that
  # spans multiple bed windows
  dup.df <- as.data.frame(olaps[dup.idx,])
  dup.spl <- split(dup.df, dup.df$subjectHits)
  dup.em <- lapply(dup.spl, function(i) {
    colMeans(as.matrix(elementMetadata(cnv[i$queryHits,])))
  })
  dup.em <- do.call(rbind, dup.em)
  
  # Initialize a metadata matrix and populate it for the BED GRanges object
  em  <- matrix(nrow=length(bed), 
                ncol=ncol(elementMetadata(cnv)), 
                dimnames = list(NULL,colnames(elementMetadata(cnv))))
  dedup.olaps <- olaps[-dup.idx,]
  em[subjectHits(dedup.olaps),] <- as.matrix(elementMetadata(cnv)[queryHits(dedup.olaps),])
  em[unique(subjectHits(olaps)[dup.idx]),] <- dup.em 
  
  # Append metadata and return
  em <- as.data.frame(em)
  bed$ID <- em$ID <- paste0("bin_", c(1:nrow(em)))
  
  return(list(seg=bed, genes=em))
}

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

reduceEsetMats <- function(gene.lrr, cols, features='SYMBOL', ord=FALSE,
                           keys=c("ENTREZ", "SYMBOL", "ENSEMBL")){
  lapply(cols, function(each.col, features){
    print(each.col)
    m <- suppressWarnings(Reduce(f=function(x,y) merge(x,y,by=keys),
                                 lapply(gene.lrr, function(i) i[['genes']][,c(keys, each.col)])))
    if(ord) m <- m[match(gene.lrr[[1]][['genes']][,keys], m[,keys]),]
    if(any(duplicated(m[,features]))) m <- m[-which(duplicated(m[,features])),]
    if(any(is.na(m[,features]))) m <- m[-which(is.na(m[,features])),]
    rownames(m) <- m[,features]
    m <- m[,-c(1:length(keys))]
    colnames(m) <- names(gene.lrr) 
    as.matrix(m)
  }, features=features)
}

##############
#### Main ####
handle <- 'gCSI'
map.to <- 'bin' #bin, gene, tad
## To add: PGA/wGII PSet (different thresholds)

switch(handle,
       'UHN'={
         anno.name <- 'UHN'
         seg.file <- 'UHN_breast.anno.seg'
         metadata <- 'ref/UHN_meta.csv'
         seg.ids <- c('sample', 'chr')
         cols <- c('CN', 'nMajor', 'nMinor')
       }, 
       'CCLE'={
         anno.name <- 'CCLE'
         seg.file <- 'CCLE.segmaf.seg'
         metadata <- 'ref/CCLE_sample_info_file_2012-10-18.txt'
         seg.ids <- c('sample', 'Chromosome')
         cols <- c('copy.ratio', 'hscr.a1', 'hscr.a2', 'modal.a1', 'modal.a2')
       }, 
       'GDSC'={
         anno.name <- 'GDSC'
         seg.file <- 'GDSC.segmaf.seg'
         seg.ids <- c('sample', 'Chromosome')
         cols <- c('copy.ratio', 'hscr.a1', 'hscr.a2', 'modal.a1', 'modal.a2')
       },
       'gCSI'={
         anno.name <- 'gCSI'
         seg.file <- c('gCSI.seg', 'gCSI_cn.raw.tsv')
         seg.ids <- c('ID', 'chrom')
         metadata <- 'ref/genentech.snparray.manifest.csv'
         cols=c("seg.mean", "nAraw", "nBraw", "nMinor", "nMajor")
       })

seg <- read.table(seg.file[1], sep="\t", header=TRUE, 
                  stringsAsFactors = FALSE, check.names = FALSE)
if(handle=='gCSI'){
  segB <- read.table(seg.file[2], sep="\t", header=TRUE, 
                     stringsAsFactors = FALSE, check.names = FALSE)
  colnames(segB)[1:4] <- c("ID", "chrom", "loc.start", "loc.end")
  seg <- merge(seg, segB, by=c("ID", "chrom", "loc.start", "loc.end"), all.x=TRUE)
}

## Annotate the CN segments
if(map.to == 'bin') windowed.bed <- genWindowedBed(bin.size=5000)

cnseg.list <- split(seg, f=seg[,seg.ids[1]])
gene.lrr <- mclapply(cnseg.list[1:5], function(cl.i){
  uid = unique(cl.i[,seg.ids[1]])
  print(paste0(uid, " - (", 
               grep(uid, names(cnseg.list)), "/", 
               length(cnseg.list), ")"))
  
  if(any(is.na(cl.i[,seg.ids[2]]))) cl.i <- cl.i[-which(is.na(cl.i[,seg.ids[2]])),]
  cnv <- makeGRangesFromDataFrame(cl.i, keep.extra.columns = TRUE, 
                                  start.field = c('Start.bp', 'loc.start', 'startpos'), 
                                  end.field = c('End.bp', 'loc.end', 'endpos'))  
  mcols(cnv) <- round(cl.i[,cols], 3)
  suppressWarnings(seqlevelsStyle(cnv) <- 'UCSC')
  
  if(map.to=='bin'){
    # Map CNV segments to a reference bed
    cl.anno <- segmentCNVs(cnv, windowed.bed)
  } else if(map.to == 'genes'){
    # Map CNV segments to genes
    cl.anno <- suppressMessages(annotateCnvs(cnv, txdb, 
                                             anno=txdb.genes,
                                             cols=cols))
    names(cl.anno) <- c('seg', 'genes')
  }
  cl.anno
})

mats <- switch(map.to,
               "bin"={
                 reduceEsetMats(gene.lrr, cols, features='ID', 
                                keys='ID', ord=TRUE)
               },
               "gene"={
                 reduceEsetMats(gene.lrr, cols, features='SYMBOL')
               })

example.seg <- gene.lrr[[1]][['seg']]
save(mats, example.seg, 
     file=file.path('esets', "raw_mats", paste0(anno.name, ".", map.to, '.Rdata')))

## Format the phenoData
if(handle == 'gCSI'){
  meta <- read.table(metadata, sep=",", header=TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
  meta$INVENTORY_SAMPLE_NAME <- gsub("Jurkat,", "Jurkat", meta$INVENTORY_SAMPLE_NAME)
  meta <- as.data.frame(meta[match(names(cnseg.list), meta$INVENTORY_SAMPLE_NAME),])
  rownames(meta) <- as.character(meta$INVENTORY_SAMPLE_NAME)
} else if (handle == 'CCLE'){
  meta <- read.table(metadata, sep="\t", header=TRUE,
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
} else if (handle == 'UHN'){
  meta <- read.csv(metadata, header=TRUE, quote="", fill=FALSE,
                   stringsAsFactors = FALSE, check.names = FALSE)
  meta <- meta[,-c(1,2)]
  dups <- which(duplicated(meta$ID_Corr))
  
  new.ids <- meta$ID_Corr
  new.ids[dups] <- paste0(new.ids[dups], "(DUP)")
  rownames(meta) <- new.ids
  
  meta <- meta[which(rownames(meta) %in% colnames(mats[[1]])),]
  meta <- meta[colnames(mats[[1]]),]
}

## Assemble the assayData environment
eset.env <- new.env()
assign("exprs", mats[[1]], envir=eset.env)
assign("A", mats[[2]], envir=eset.env)
assign("B", mats[[3]], envir=eset.env)
if(handle != 'UHN'){
  assign("modalA", mats[[4]], envir=eset.env)
  assign("modalB", mats[[5]], envir=eset.env)
}




## Assemble the eset
cl.phenoData <- new("AnnotatedDataFrame", data=meta)
cl.eset <- ExpressionSet(assayData=eset.env,
                         phenoData=cl.phenoData,
                         annotation=anno.name)
save(cl.eset, file=file.path("esets", paste0(anno.name, "_eset.Rdata")))
