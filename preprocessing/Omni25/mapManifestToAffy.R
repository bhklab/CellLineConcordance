# Date: May-26-2018
# Author: Rene Quevedo
# Purpose: Takes the trimmed versions of Omni2.5 Quad manifest and the Affymetrix
#   SNP6.0 csv annotations and maps the Omni to the Affy array.  It will find the
#   SNPs found in both datasets by matching genomic locations, then strand-correct
#   the SNPs, followed by identifying which allele is "Allele A" and which is 
#   "Allele B".  Results are saved in an Annotated Dataframe

library(magrittr)
library(Biobase)
#' getSnp
#'
#' @param x Vector of GenomeStudio format SNPs, (e.g. [A/T], [C/T], [T/A])
#' @param snp Either "A", or "B" to return either or
#'
#' @return If A, the return will be A, C, T.  Likewise, if B, the return 
#' will be T, T, A
#' @export
#'
#' @examples
getSnp <- function(x, snp){
  if(snp=='A'){
    x %<>%
      gsub('\\[', "", .) %>%
      gsub("/.*", "", .) 
  } else if (snp =='B'){
    x %<>%
      gsub('\\[./', "", .) %>%
      gsub("\\]", "", .) 
  } else {
    stop("Specifiy snp=A or snp=B")
  }
  x
}

#' syncStrandedSnps
#'
#' @param snp Vector of alleles (e.g. A, A, T, C)
#' @param strand_x Vector of reference strand (e.g. +, -, +, -)
#' @param strand_y Vector of comparing strand (e.g. -, -, -, -)
#'
#' @return Vector of SNPs that are reverse complemented when the comparing strand
#' doesnt match the reference strand
#' @export
#'
#' @examples
syncStrandedSnps <- function(snp, strand_x, strand_y){
  revcomp <- list("A"="T",
                  "T"="A",
                  "C"="G",
                  "G"="C")
  strand.idx <- which(strand_x != strand_y)
  snp[strand.idx] <- sapply(snp[strand.idx], function(x) revcomp[[x]])
  snp
}

#' setSnpToggleFlag
#'
#' @param snp_x1 Vector for Reference dataset Allele A (e.g. C, C, T)
#' @param snp_x2 Vector for Reference dataset Allele B (e.g. A, T, G)
#' @param snp_y1 Vector for Comparing dataset Allele A (e.g. A, T, T)
#' @param snpy_y2 Vector for Comparing dataset Allele B (e.g. C, C, G)
#'
#' @return A Flag of NA, CORRECT or FLIP to indicate whether A-B of reference
#' matches A-B of comparing dataset (CORRECT) or B-A of comparing dataset (FLIP)
#' @export
#'
#' @examples
setSnpToggleFlag <- function(snp_x1, snp_x2, snp_y1, snp_y2){
  tot.orient <- rep(NA, length(snp_x1))
  
  corr.orient <- (snp_x1 == snp_y1) & (snp_x2 == snp_y2)
  tot.orient[which(corr.orient)] <- 'CORRECT'
  
  flip.orient <- (snp_x1 == snp_y2) & (snp_x2 == snp_y1)
  tot.orient[which(flip.orient)] <- 'FLIP'
  
  tot.orient
}


############
##  MAIN
############
setwd("/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/reference/annotations/mapping")
affy <- read.table("GenomeWideSNP_6.na35.annot.trimmed.csv", header=TRUE, stringsAsFactors = FALSE, check.names=FALSE)
omni <- read.table("HumanOmni2.5-4v1_H.trimmed.csv", header=TRUE, stringsAsFactors = FALSE, check.names=FALSE)

affy.id <- apply(affy[,c("Chromosome", "Physical_Position")], 1, function(x)
  paste(x, collapse=","))
omni.id <- apply(omni[,c("Chr", "MapInfo")], 1, function(x)
  paste(x, collapse=","))
omni.id <- gsub(" ", "", omni.id)

affy$mergeid <- affy.id
omni$mergeid <- omni.id

affy.omni <- merge(x=affy, y=omni, by="mergeid", all=TRUE)
affy.omni.full <- affy.omni[which(apply(affy.omni, 1, function(x) !any(is.na(x)))),]

affy.omni.full$SNP_A <- getSnp(affy.omni.full$SNP, 'A')
affy.omni.full$SNP_B <- getSnp(affy.omni.full$SNP, 'B')

affy.omni.full$scSNP_A <- with(affy.omni.full, syncStrandedSnps(SNP_A, Strand, RefStrand))
affy.omni.full$scSNP_B <- with(affy.omni.full, syncStrandedSnps(SNP_B, Strand, RefStrand))

affy.omni.full$flip <- with(affy.omni.full, setSnpToggleFlag(Allele_A, Allele_B, 
                                                             scSNP_A, scSNP_B))
metaData <- data.frame(labelDescription=c(
  "Unique merge ID ([chr],[pos])",
  "Affy6 probe set IDs",
  "Affy6 dbSNP IDs",
  "Affy6 Chromosome",
  "Affy6 Genomic Loci (GRCh37)",
  "Affy6 Strand", 
  "Affy6 A-allele",
  "Affy6 B-allele",
  "Omni2.5quad probe set IDs",
  "Omni2.5quad SNPs",
  "Omni2.5quad Chromosome",
  "Omni2.5quad Genomic Loci (GRCh37)",
  "Omni2.5quad Strand",
  "Omni2.5quad A-allele",
  "Omni2.5quad B-allele",
  "Omni2.5quad strand corrected A-allele (ref=Affy6)",
  "Omni2.5quad strand corrected b-allele (ref=Affy6)",
  "Omni2.5 order of alleles in reference to Affy6"))
affy.omni.full <- AnnotatedDataFrame(data=affy.omni.full, varMetadata=metaData)

save(affy.omni.full, file="affy_omni.RData")
write.table(pData(affy.omni.full), file="affy_omni.tsv",
            col.names = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)
