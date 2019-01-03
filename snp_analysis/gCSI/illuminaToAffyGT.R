###############################################################
#
#                 Illumina OMNI2.5 To Affy6 Genotype Converter
#
# Author: Rene Quevedo
# Date Created: May 31-2016
###############################################################
# Function: Takes a CRLMM genotype datastructure and converts to
#   the genotype format 0/0,
#   0/1, and 1/1 to the birdseed genotype standards: 0, 1, 2.
#   This is done for all SNPs that overlap with the reference 
#   Affymetrix SNP 6 annotation file
# VCF Files Tested:
#   -Illumina HumanOmni 2.5-8 Quad - GNE Dataset
#
#   Usage: Rscript illuminaToAffyGT.R [illuminaTsv] [illuminaAnno] [affyAnno] [outfile]
###############################################################
args = commandArgs(trailingOnly=TRUE)
tsvFile <- args[1]
affy6Anno <- args[2]
omniAnno <- args[3]
outFile <- args[4]



#######  FUNCTIONS
show.functions <- 1
if(show.functions == 1){
  #Taking into consider Strand direction, returns the Birdseed Genotypes for any given
  #VCF genotype
  #   For example: Convert Chr1:1234C/T  vcf genotype: 1/1 to affy6 genotype 2
  getVcfGeno <- function(x){
    #Get VCF Genotype and validate:
    genotype <- gsub(":.+", "", x[10])
    affy.genotype <- -1
    if(length(grep("[^10/]", genotype)) > 0){
      print(paste("Error: Chr", x['CHROM'], " - ", x['POS'], " - Strand", x['Strand'], sep=""))
      print(paste("Script currently does not account for complex genotypes: ", genotype, sep=""))
    } else {
      affy.genotype <- convertVcfToAffyGT(genotype)
    }
    
    # If strand is reverse, return the complement allele
    complementAllele <- c(A="T", T="A", C="G", G="C")
    if(x['Strand'] %in% '-'){
      x['Allele_A'] <- complementAllele[x['Allele_B']]
      x['Allele_B'] <- complementAllele[x['Allele_A']]
    }
    
    # Ensures that the Ref allele in our vcf matches Allele_A in affy6
    complementGenotype <- c('0'='2', '1'='1', '2'='0')
    allele.order <- match(x['Allele_A'], c(x['REF'], x['ALT']))
    if (is.na(allele.order)){
      print(paste("Error: Chr", x['CHROM'], " - ", x['POS'], " - Strand", x['Strand'], sep=""))
      print(paste("Could not match allele order between ",
                  x['REF'], "/", x['ALT'], " and ",
                  x['Allele_A'], "/", x['Allele_B'], sep=""))
      affy.genotype <- -1
    } else if (allele.order == 2){
      ref.allele <- x['REF']
      x['REF'] <- x['ALT']
      x['ALT'] <- ref.allele
      
      affy.genotype <- complementGenotype[as.character(affy.genotype)]
      affy.genotype <- as.integer(affy.genotype)
    } 
    
    return(affy.genotype)
  }
  
  # Converts VCF Genotype to Affy6 Genotype
  convertVcfToAffyGT <- function(vcf.gt){
    affy.gt <- -1
    if(vcf.gt %in% '1/1'){
      affy.gt <- 2
    } else if((vcf.gt %in% '0/1') || (vcf.gt %in% '1/0')){
      affy.gt <- 1
    } else if(vcf.gt %in% '0/0'){
      affy.gt <- 0
    } 
    return(affy.gt)
  }
}

#######  MAIN
# Read in Data and format:
vcf.df <- read.csv(tsvFile, sep="\t", header=TRUE, comment.char = "#",
                   stringsAsFactors=FALSE, check.names=FALSE)
vcf.df$CHROM <- gsub("^chr", "", vcf.df$CHROM, ignore.case = TRUE)
vcf.pos <- paste(vcf.df$CHROM, vcf.df$POS, sep=",")

affyAnno.df <- read.csv(affy6Anno, sep="\t", header=TRUE, comment.char = "#",
                        stringsAsFactors=FALSE, check.names=FALSE)
affyAnno.df$Chromosome <- gsub("^chr", "", affyAnno.df$Chromosome, ignore.case=TRUE)
affyAnno.pos <- paste(affyAnno.df$Chromosome, affyAnno.df$Physical_Position, sep=",")

# Find the Affy6 Probe Id for matching loci:
affyMatch.idx <- match(vcf.pos, affyAnno.pos)
vcf.idx <- which(!is.na(affyMatch.idx))

vcf.anno <- affyAnno.df[affyMatch.idx,]
vcf.pos.anno <- cbind(vcf.df[vcf.idx, ], vcf.anno[vcf.idx,])

#Align the REF and ALT allele to match for the Affy6 Allele_A and Allele_B (factor in strand direction)
vcf.affyGT <- apply(vcf.pos.anno, 1, getVcfGeno)
vcf.pos.anno$vcf.affyGT <- vcf.affyGT
vcf.pos.anno <- vcf.pos.anno[which(vcf.pos.anno$vcf.affyGT != -1),]

# Writes the Genotype to file for further matching:
vcf.affy.df <- vcf.pos.anno[c('Probe_Set_ID', 'vcf.affyGT')]
colnames(vcf.affy.df) <- c('probeset_id', gsub(".vcf", "", basename(tsvFile)))
write.table(vcf.affy.df, file = outFile, quote=FALSE, sep="\t", 
            row.names=FALSE, col.names=TRUE)

