# Date: May-26-2018
# Author: Rene Quevedo
# Purpose: Syncs up the GNE manifest obtained from the Genetech people and the Unknown IDs 
#   that were assigned during run of GenomeStudio
library(magrittr)

gne_manifest <- '/mnt/work1/users/bhklab/Data/GNE/genentech.snparray.manifest.csv'
gne_samplesheet <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2017_GNE/illumina/Ref/samplesheet.tmp'

gne_manifest <- read.table(gne_manifest, header = TRUE, sep=",", stringsAsFactors = FALSE, check.names = FALSE)
gne_samplesheet <- read.table(gne_samplesheet, header = TRUE, sep=",", stringsAsFactors = FALSE, check.names = FALSE)

gne_samplesheet$IDATid <- apply(gne_samplesheet, 1, function(x) paste(c(x['SentrixBarcode_A'], 
                                                                        x['SentrixPosition_A']), 
                                                                      collapse="_"))

corrRandC <- function(x){
  x %<>%
    gsub("r", "R", .) %>%
    gsub("c", "C", .)
  x
}
gne_samplesheet$IDATid <- corrRandC(gne_samplesheet$IDATid)
gne_manifest$IDATid <- corrRandC(gne_manifest$IDATid)

gne_mapping <- merge(gne_samplesheet, gne_manifest, 
                     by.x="IDATid", by.y="IDATid", all=TRUE)
na.idx <- is.na(gne_mapping$Sample_ID)
gne_mapping.na <- gne_mapping[which(na.idx),]
gne_mapping <- gne_mapping[which(!na.idx),]
gne_mapping <- rbind(gne_mapping, gne_mapping.na)
write.table(gne_mapping, file="GNE_mapping.tsv", sep="\t", 
            col.names = TRUE, quote = FALSE, row.names = FALSE)