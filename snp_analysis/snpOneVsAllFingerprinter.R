.libPaths(c(.libPaths(), "/mnt/work1/users/pughlab/bin/r_lib/library/3.1"))
###############################################################
#
#                 SNP Fingerprinter One-vs-All
#
# Author: Rene Quevedo
# Date Created: August 25-2016
###############################################################
# Function: Creates and plots a one-against-all matrix of
#   sample vs cell lines from all datasets .  
#   Creates a heatmap of percent identity.
###############################################################
library(gplots)

###############################
#         Variables
###############################
args <- commandArgs(trailingOnly = TRUE)

refGeno <- args[1]
sampleGeno <- args[2]
cclAnno <- args[3]
outFile <- args[4]


###############################
#         Functions
###############################
# Function: calcPercIdentity
# Purpose:  Given two sets of ordered genotypes from birdseed, it will
#   calculate the percentage identity between the two genotypes.
# Input:  cl.call.a <- cell-line genotype calls A
#         cl.call.b <- cell-line genotype calls B
#         cl.conf.a <- cell-line genotype conf A
#         cl.conf.b <- cell-line.genotype conf B
# Returns:  String Vector: Statement of how to use the script
calcPercIdentity <- function(cl.call.a, cl.call.b,
                             cl.conf.a = NULL, cl.conf.b = NULL){
  num.of.raw.probes <- length(cl.call.a)
  num.of.matches <- length(cl.call.a[cl.call.a == cl.call.b])  
  perc.identity <- num.of.matches / num.of.raw.probes 
  
  return(perc.identity)
}


###############################
#           Main
###############################

#Load the dataframes
ref.df <- read.table(refGeno, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
sample.df <- read.table(sampleGeno, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
load(cclAnno)

# Creates a matrix of all genotype identity matches
sample.df <- sample.df[match(ref.df$probeset_id, sample.df$probeset_id),]
ref.df <- ref.df[,-which(colnames(ref.df) %in% 'probeset_id')]
genotype.conc <- apply(ref.df, 2, function(j) calcPercIdentity(sample.df[[2]],j))
genotype.conc <- as.data.frame(genotype.conc)
colnames(genotype.conc) <- "concordance"

cl.id <- list()
dataset.ids <- sapply(rownames(genotype.conc), function(x) grep(x, cell.line.anno, ignore.case=TRUE))
dataset.ids <- as.integer(sapply(dataset.ids, function(x) x[1]))
for(each.ds.id in c(1:length(dataset.ids))){
  if(!is.na(dataset.ids[each.ds.id])){
    row.idx <- grep(rownames(genotype.conc)[each.ds.id], 
                    cell.line.anno[,dataset.ids[each.ds.id]], 
                    ignore.case=TRUE)
    data.set <- colnames(cell.line.anno)[dataset.ids[each.ds.id]][1]
    cl <- cell.line.anno[row.idx, 'unique.cellid'][1]
  } else {
    data.set <- NA
    cl <- rownames(genotype.conc)[each.ds.id]
  }
  
  cl.df <- data.frame(dataset=data.set, cellId=cl, sampleId=basename(sampleGeno))
  cl.id[[each.ds.id]] <- cl.df
}

cl.id.df <- do.call("rbind", cl.id)
genotype.conc <- cbind(genotype.conc, cl.id.df)

write.table(genotype.conc, file=outFile, sep="\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
print("Completed genotype identity matching. Outputting Rdata file..")





