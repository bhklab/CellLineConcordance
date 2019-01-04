# Date: May-26-2018
# Author: Rene Quevedo
# Purpose: 
library(magrittr)
library(Biobase)

#' annotateConc
#'
#' @param omni.conc Ref Samples x Omni sample concordance matrix
#'
#' @return The top matches and annotated matches for each Omni sample.
#'  A 2nd column of 1 and 0 describes "anno-match" and "threshold-match" respectively
#' @export
#'
#' @examples
annotateConc <- function(omni.conc, threshold=0.80){
  lapply(seq_along(colnames(omni.conc)), function(x){
    match.idx <- grep(paste0("_", colnames(omni.conc)[x], "$"), rownames(omni.conc))
    threshold.idx <- which(omni.conc[,x] >= threshold)
    all.idx <- unique(sort(c(match.idx, threshold.idx)))
    match.thresh <- rep(0, length(all.idx))
    match.thresh[!is.na(match(all.idx, match.idx))] <- 1
    
    all.matches <- omni.conc[all.idx, x, drop=FALSE]
    all.matches <- cbind(all.matches, match.thresh)
    all.matches[order(all.matches[,2], decreasing = TRUE), ,drop=FALSE]
  })
}

#' flipAlleles
#'
#' @param omni.tmp Allele matrix (omni probe x samples) containing only AA, AB, or BB
#' @param snp.mapping The mapping of Omni snps to Affy snps with only FLIP or CORRECT in flip column
#'
#' @return The flipped omni alleles matrix
#' @export
#'
#' @examples
flipAlleles <- function(omni.tmp, snp.mapping){
  # Creates a boolean matrix indicating which rows to flip the alleles
  boolean.flip.mat <- matrix(FALSE, ncol=ncol(omni.tmp), nrow=nrow(omni.tmp))
  flip.boolean <- snp.mapping$flip == 'FLIP'
  boolean.flip.mat[flip.boolean,] <- TRUE
  
  # Uses gsub to flip "BB"->"AA" and "AA"->"BB"
  omni.tmp[boolean.flip.mat] %<>%
    gsub("BB", "0", .) %>%
    gsub("AA", "2", .) %>%
    gsub("0", "AA", .) %>%
    gsub("2", "BB", .) 
  
  omni.tmp
}

#' genotypeConcordance
#'
#' @param omni.tmp A Matrix of SNPs by Samples for omni data in 0, 1, 2 format
#' @param ref.tmp A Matrix of SNPs by Samples for affy data in 0, 1, 2 format
#'
#' @return A matrix of concordances, Affy x Omni samples, in a range of 0 to 1
#' @export
#'
#' @examples
genotypeConcordance <- function(omni.tmp, ref.tmp){
  omni.conc <- apply(omni.tmp, 2, function(x){
    print(head(x))
    x <- as.integer(as.character(x))
    x <- matrix(rep(x, ncol(ref.tmp)), ncol=ncol(ref.tmp))
    xy <- apply(ref.tmp == x, 2, sum)
    round(xy / nrow(ref.tmp),3)
  })
  
  omni.conc
}

#' mapIds
#'
#' @param cell.line.anno Cell line annotation dataframe
#' @param id.mapping Mapping between GenomeStudio "Unk" ids and the actuall cell IDs
#' @param omni.conc Omni x Affy sample concordance matrix between 0-1
#' @param ref.tmp SNP x Affy sample allele matrix with colnames set as .CEL files
#' @param set.type Either 'rownames' or 'colnames' to set the corresponding value for omni.conc
#'
#' @return A vector of cell line IDs
#' @export
#'
#' @examples
mapIds <- function(cell.line.anno, id.mapping, omni.conc, ref.tmp=NULL, set.type){
  if(set.type=='colnames'){
    unk.ids <- gsub(".GType|.n[AB]raw.Rdata", "", colnames(omni.conc))
    map.ids <- id.mapping[match(unk.ids, id.mapping$Sample_ID), 'INVENTORY_SAMPLE_NAME']
  } else if (set.type=='rownames'){
    if(!is.null(ref.tmp)){
      ref.ids <- colnames(ref.tmp)
    } else {
      ref.ids <- rownames(omni.conc)
    }
    
    map.ids <- sapply(ref.ids, function(x) {
      match.idx <- which(x == cell.line.anno, arr.ind=TRUE)
      if(nrow(match.idx)>1) match.idx <- match.idx[1, ,drop=FALSE]
      cell.id <- cell.line.anno[match.idx[,1],'unique.cellid']
      paste0(gsub(".filename.*", "", colnames(cell.line.anno)[match.idx[,2]]), 
             "_", cell.id)
    })
  }
  map.ids
}
###############
##  MAIN
###############
setwd("/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2017_GNE/snp_fingerprint")
args = commandArgs(trailingOnly=TRUE)
ref.mat <- file.path("input", "ccle.cgp.gdsc.pfizer.birdseed-v2.calls.txt")
omni.mat <- file.path("input", paste0("GNE_Genotype_", args[1], ".txt"))
snp.mapping <- file.path("reference", "affy_omni.RData")
id.mapping <- file.path("reference", "GNE_mapping.tsv")
anno.mapping <- file.path("reference", "merged_annotations.Rdata")
#save(ref.mat, file="/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/rdataFiles/snp/ccle.cgp.gdsc.pfizer.birdseed-v2.calls.RData")
#save(ref.mat, file="/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/rdataFiles/snp/ccle.cgp.gdsc.pfizer.birdseed-v2.omniOverlap.RData")

#ref.mat <- read.table(ref.mat, header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
#load("/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/rdataFiles/snp/ccle.cgp.gdsc.pfizer.birdseed-v2.calls.RData")
load("/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/rdataFiles/snp/ccle.cgp.gdsc.pfizer.birdseed-v2.omniOverlap.RData")
load(anno.mapping)
omni.mat <- read.table(omni.mat, header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
load(snp.mapping)
snp.mapping <- pData(affy.omni.full); rm(affy.omni.full)
snp.mapping <- snp.mapping[which(!is.na(snp.mapping$flip)), ]
id.mapping <- read.table(id.mapping, header=TRUE, sep=",", stringsAsFactors = FALSE, check.names = FALSE)

# Prints Dimensions
print(dim(omni.mat))
print(dim(ref.mat))
print(dim(snp.mapping))

# Subsets the Omni to Affy6 SNP.mapping dataframe to only those 
# that are found in the Affy6 SNP set. Reorders to match the Affy6 order as well.
affy.probes <- rownames(ref.mat)
omni.probes <- rownames(omni.mat)
affy.keep.idx <- which(affy.probes %in% snp.mapping$Probe_Set_ID)
included.probes.idx <- snp.mapping$Probe_Set_ID %in% affy.probes[affy.keep.idx]
snp.mapping <- snp.mapping[which(included.probes.idx),]
snp.mapping <- snp.mapping[match(affy.probes[affy.keep.idx], snp.mapping$Probe_Set_ID),]

# Reorders the Omni SNP matrix to match that of Affy6 SNP and SNP mapping order
if(length(affy.keep.idx) != nrow(ref.mat)) ref.mat <- ref.mat[affy.keep.idx, -1] 
print(head(match(snp.mapping$Name, omni.probes)))
omni.mat <- omni.mat[match(snp.mapping$Name, omni.probes), -c(1,2)]

if(all(rownames(ref.mat) == snp.mapping$Probe_Set_ID)){
  print("All affy are ordered...")
  if(all(rownames(omni.mat) == snp.mapping$Name)){
    print("All omni are ordered...")
  } else {
    tf <- rownames(omni.mat) == snp.mapping$Name
    print(table(tf))
    print(head(rownames(omni.mat)[which(!tf)]))
    print(head(snp.mapping$Name[which(!tf)]))
    print("Omni not ordered.")
  }
} else {
  "Affy not ordered."
}
print(dim(snp.mapping))
print(dim(omni.mat))

omni.mat <- flipAlleles(omni.mat, snp.mapping)
omni.mat[matrix(TRUE, ncol=ncol(omni.mat), nrow = nrow(omni.mat))] %<>%
  gsub("AA", "0", .) %>%
  gsub("AB", "1", .) %>%
  gsub("BB", "2", .) %>%
  gsub("NC", "-1", .)


omni.conc <- genotypeConcordance(omni.mat[,1:5], ref.mat)
colnames(omni.conc) <- mapIds(cell.line.anno, id.mapping, omni.conc, ref.mat, 'colnames')
rownames(omni.conc) <- mapIds(cell.line.anno, id.mapping, omni.conc, ref.mat, 'rownames')


top.conc <- annotateConc(omni.conc)
names(top.conc) <- colnames(omni.conc)
save(omni.conc, top.conc, file=file.path("output", paste0("GNE_Conc_", args[1], ".Rdata")))
  
