#' getGenes
#'
#' @param genome.build ['character'] either hg18, hg19, or hg38
#'
#' @return A GRanges object id with ENTREZ IDs for genes
#' @export
#'
#' @examples
getGenes <- function(genome.build="hg19"){
  switch(genome.build,
         hg18={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg18.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg18.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg18.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg18.knownGene
         },
         hg19={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg19.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg19.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg19.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg19.knownGene
         },
         hg38={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg38.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg38.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg38.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg38.knownGene
         },
         stop("genome must be 'hg19' or 'hg38'"))
  
  genes0 <- genes(package)
  idx <- rep(seq_along(genes0), elementNROWS(genes0$gene_id))
  genes <- granges(genes0)[idx]
  genes$gene_id = unlist(genes0$gene_id)
  genes
}

#' annotateSegments
#'
#' @param cn.data Copy-number data as a data.frame
#' @param genes A GRanges object from getGenes() function
#' @param mart A biomaRt object
#' @param use.mart [boolean] whether to load a mart if none is given
#'
#' @return a list containing a GRanges object with annotated HUGO genes
#'  and a dataframe containing the HUGO gene and their seg.mean
#' @export
#'
#' @examples
annotateSegments <- function(cn.data, genes, mart=NULL, use.mart=FALSE){
  suppressPackageStartupMessages(require(biomaRt))
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  suppressPackageStartupMessages(require(GenomicRanges))
  suppressPackageStartupMessages(require(AnnotationDbi))
  if(use.mart & is.null(mart)){
    mart <- useMart("ENSEMBL_MART_ENSEMBL")
    mart <- useDataset("hsapiens_gene_ensembl", mart)
  }
  
  
  gr0 <- makeGRangesFromDataFrame(cn.data,keep.extra.columns=TRUE)
  seqlevelsStyle(gr0) <- 'UCSC'
  olaps <- findOverlaps(genes, gr0, type="within")
  idx <- factor(subjectHits(olaps), levels=seq_len(subjectLength(olaps)))
  gr0$gene_ids <- splitAsList(genes$gene_id[queryHits(olaps)], idx)
  gr0$gene_ids <- lapply(gr0$gene_ids, function(input.id) {
    if(length(input.id) > 0){ 
      tryCatch({
        ens <- mapIds(org.Hs.eg.db,
                      keys=input.id,
                      column="SYMBOL",
                      keytype="ENTREZID",
                      multiVals="first")
        
        if(!is.null(mart)){
          ens[which(is.na(ens))] <- getBM(
            mart=mart,
            attributes="external_gene_name",
            filter="entrezgene",
            values=names(ens[which(is.na(ens))]),
            uniqueRows=TRUE)
        }
        ens
      }, error=function(e){NULL})
    } else { NA }
  })
  return(gr0)
}
