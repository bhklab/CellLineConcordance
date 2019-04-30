library(weights)
library(CollapsABEL)
library(grid)
library(optparse)
library(bigmemory)
library(biganalytics)

####################
##### SET-UP
option_list = list(
  ## Reference CNV/StdDev Paths [r,m,s,t]
  make_option(c("-r", "--refcnvdir"), type="character",
                default='/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/rdataFiles/cnv',
              help="Directory containing the CNV matrices for CCLE/GDSC/CGP/Pfizer [default=%default]",
              metavar="character"),
  make_option(c("-m", "--metadataref"), type="character",
              default='/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/reference',
              help="Directory containing cell line metadata for the reference datasets [default=%default]",
              metavar="character"),
  make_option(c("-s", "--refcndir"), type="character",
              default='v2_CNsegments',
              help="Directory containing the CNV bigmemory matrices for Ref datsets [default=%default]",
              metavar="character"),
  make_option(c("-t", "--refvardir"), type="character",
              default='v1_20xStdev',
              help="Directory containing the StdDev list for Ref datsets [default=%default]",
              metavar="character"),
  
  ## Cytoband and SNP dir [b, n, c]
  make_option(c("-b", "--cytodir"), type="character",
              default='~/git/reference/cytoband',
              help="Directory containing the UCSC hg19 cytoband reference files [default=%default]",
              metavar="character"),
  make_option(c("-n", "--snpdir"), type="character",
              default='~/git/snp_fingerprint/environments',
              help="Directory containing the SNP concordance matrix [default=%default]",
              metavar="character"),
  make_option(c("-c", "--snpmatrix"), type="character",
              default='plotPairwiseSnp.v2.Rdata',
              help="SNP concordance matrix file name [default=%default]",
              metavar="character"),
  
  ## Alternate CNV/StdDev Paths [x,o,a]
  make_option(c("-x", "--altcnvdir"), type="character",
              default=NULL,
              help="Directory containing the CNV matrices for the comparing CNV matrix (e.g. GNE) [default=%default]",
              metavar="character"),
  make_option(c("-o", "--outdir"), type="character",
              default="/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/cnv_analysis/output",
              help="Output directory [default=%default]",
              metavar="character"),
  make_option(c("-a", "--annoref"), type="character",
              default=NULL,
              help="Ref datasets annotation metadata [default=%default]",
              metavar="character"),
  
  ## Functions Path
  make_option(c("-f", "--functions"), type="character",
              default='~/git/CellLineConcordance/cnv_analysis/src',
              help="Source directory containing all functions used [default=%default]",
              metavar="character"),
  make_option(c("-v", "--vmem"), type="character",
              default=TRUE,
              help="Uses a local bigmemory object for vmem [default=%default]",
              metavar="logical")
)
opt_parser = OptionParser(usage = "Rscript %prog [options]", option_list=option_list)
opt=parse_args(opt_parser)

###############################
#     Interactive Workspace
interactive=TRUE
if (interactive) {
  opt$altcnvdir <- "/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/rdataFiles/cnv/GNE"
  opt$annoref <- "/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2017_GNE/snp_fingerprint/reference/GNE_mapping.tsv"
}

###############################
#     Environment
## Reference
ccl.cnv.datadir <- opt$refcnvdir
seg.dir <- opt$refcndir
stdev.dir <- opt$refvardir
out.dir <- opt$outdir
anno.dir <- opt$metadataref
anno.file <- "merged_annotations.Rdata"

cyto.dir <- opt$cytodir
snp.matrix.dir  <- opt$snpdir
snp.matrix.file <- opt$snpmatrix
cyto.file="chr_cytoband.hg19.Rdata"

## Alternate
id.mapping <- opt$annoref
alt.seg.dir <- opt$altcnvdir

load(file.path(cyto.dir, cyto.file))  # "cytoband.df", "raw.cytoband.df", "chr.df"  
load(file.path(snp.matrix.dir, snp.matrix.file))  # "y.df", "low.match.list", "cell.line.anno", "matrix.perc.match"
load(file.path(anno.dir,anno.file))


setwd(out.dir)
threshold <- 0.8
max.length <- 10
hscr <- 'hscra2'

fn.headers <- c('cgp.filename', 'ccle.filename', "pfizer.filename.x", "pfizer2.filename.y", 'pfizer3.filename', 'gdsc.filename.x', 'gdsc2.filename.y')
names(fn.headers) <- c("cgp", "ccle", "pfizer", "pfizer2", "pfizer3", "gdsc", "gdsc2")

ds.col <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")
names(ds.col) <- c("GDSC", "CGP", "CCLE", "PFIZER", "GNE")

###############################
#     Function
source(file.path(opt$functions, "cnvPlotFunctions.R"))


###############################
#     Main: Loading in CN-profiles
#Load in all data files
load(paste(hscr, "combined_match_mismatch.Rdata", sep="."))
mm.list <- combined.match.mismatch.list
rm(combined.match.mismatch.list)


## Load in Reference CNVs/StdDevs Files
#writeBigMatrix(hscr.mat.a1, uniq.id="hscra1", profile.id="cnProfile")
#writeBigMatrix(hscr.mat.a2, uniq.id="hscra2", profile.id="cnProfile")
#hscr.mat.a12 <- hscr.mat.a1 + hscr.mat.a2
#writeBigMatrix(hscr.mat.a12, uniq.id="hscr_a12", profile.id="cnProfile")
#stdev.mat.a12 <- do.call("cbind", stdev.list.a12)
#writeBigMatrix(stdev.mat.a12, uniq.id="hscr_a12", profile.id="StdDev")

stdev.list.a1 <- loadStdev("hscra1")
stdev.list.a2 <- loadStdev("hscra2")
stdev.list.a12 <- lapply(names(stdev.list.a1), function(x) {
  col.id.a2 <- colnames(stdev.list.a2[[x]])
  col.id.a1 <- colnames(stdev.list.a1[[x]])
  cola1.idx <- match(col.id.a2, col.id.a1)
  x <- stdev.list.a1[[x]][,cola1.idx] + stdev.list.a2[[x]]
  x / 2
})
names(stdev.list.a12) <- names(stdev.list.a1)

if(opt$vmem){
  print("Loading matrices into vmem")
  hscra1.list <- loadHscrMat(file.path(ccl.cnv.datadir, seg.dir), 'hscra1', 'mem.cnProfile')
  hscr.mat.a1 <- hscra1.list[['hscr']]
  ord.df <- hscra1.list[['ord.df']]
  hscr.mat.a2 <- loadHscrMat(file.path(ccl.cnv.datadir, seg.dir), 'hscra2', 'mem.cnProfile')
  hscr.mat.a12 <- loadHscrMat(file.path(ccl.cnv.datadir, seg.dir), 'hscr_a12','cnProfile')
  stdev.mat.a12 <- loadHscrMat(file.path(ccl.cnv.datadir, seg.dir), 'hscr_a12','StdDev')
} else {
  print("Loading matrices into memory")
  # saveMatrix <- function(bm.mat, path, file.id){
  #   hscr.mat <- round(as.matrix(hscr.mat.a1),2)
  #   save(hscr.mat, file=file.path(path, file.id))
  # }
  # saveMatrix(hscr.mat.a1, file.path(ccl.cnv.datadir, seg.dir, "round"), "hscra1.mem.cnProfile.Rdata")
  # saveMatrix(hscr.mat.a2, file.path(ccl.cnv.datadir, seg.dir, "round"), "hscra2.mem.cnProfile.Rdata")
  # saveMatrix(hscr.mat.a12, file.path(ccl.cnv.datadir, seg.dir, "round"), "hscra12.mem.cnProfile.Rdata")
  # save(ord.df, file=file.path(ccl.cnv.datadir, seg.dir, "round", "ord.df.Rdata"))
  
  load(file.path(ccl.cnv.datadir, seg.dir, "hscra1.mem.cnProfile.Rdata"))
  load(file.path(ccl.cnv.datadir, seg.dir, "hscra2.mem.cnProfile.Rdata"))
  load(file.path(ccl.cnv.datadir, seg.dir, "hscra12.mem.cnProfile.Rdata"))
  load(file.path(ccl.cnv.datadir, seg.dir, "ord.df.Rdata"))
  stdev.mat.a12 <- do.call("cbind", stdev.list.a12)
}


###############################
#     Main: Generating Results

## Load in Comparing CNVs matrices for GNE (Fig.6)
setwd(out.dir)
if(!is.null(alt.seg.dir)){
  id.mapping <- read.table(id.mapping, header=TRUE, sep=",", stringsAsFactors = FALSE, check.names = FALSE)
  
  if(opt$vmem){
    x.hscr.mat.a2 <- loadHscrMat(alt.seg.dir, hscr="nAraw", profile.id="cnProfile", id.mapping)
    x.hscr.mat.a1 <- loadHscrMat(alt.seg.dir, hscr="nBraw", profile.id="cnProfile", id.mapping)
    x.hscr.mat.a12 <- loadHscrMat(alt.seg.dir, hscr="copyratio", profile.id="cnProfile", id.mapping)
  } else {
    # saveMatrix(x.hscr.mat.a1, alt.seg.dir, "nAraw.cnProfile.Rdata")
    # saveMatrix(x.hscr.mat.a2, alt.seg.dir, "nBraw.cnProfile.Rdata")
    # saveMatrix(x.hscr.mat.a12, alt.seg.dir, "nABraw.cnProfile.Rdata")
    load(file.path(alt.seg.dir, "nAraw.cnProfile.Rdata"))
    load(file.path(alt.seg.dir, "nBraw.cnProfile.Rdata"))
    load(file.path(alt.seg.dir, "nABraw.cnProfile.Rdata"))
  }
  list.x <- list(list("i"='PSN1', "c1"='CCLE',
                      "j"='PSN1', "c2"='GNE'),
                 list("i"='NCI-H2052', "c1"='CCLE',
                      "j"='NCI-H2052', "c2"='GNE'),
                 list("i"='MDA-MB-231', "c1"="CCLE",
                      "j"='MDA-MB-231',"c2"="GNE"),
                 list("i"='SW403', "c1"='CGP',
                      "j"='SW403', "c2"='GNE'),
                 list("i"='SCC-15', "c1"='CCLE',
                      "j"='SCC-15', "c2"='GNE'),
                 list("i"='PSN1', "c1"='GDSC',
                      "j"='PSN1', "c2"='GNE'))
                 #list("i"='MDA−MB−468', "c1"='GDSC',
                #      "j"='MDA-MB-468 (IVCC)', "c2"='GNE'))
  setwd(out.dir)
  diff.x <- runAndPlotSummary(list.x, "gne", analysis.status='between')
}


## Load in Comparing CNVs matrices for Reference Datasets (Fig. 3)
genfig3 <- TRUE
if(genfig3){
  discordant.list <- list(list("i"='HPAC', "c1"="CCLE",
                               "j"='KCI-MOH1',"c2"="CCLE"),
                          list("i"='OCI-LY10', "c1"="CCLE",
                               "j"='KARPAS-422',"c2"="CCLE"))
  drift.list <- list(list("i"='NCI-H23', "c1"="GDSC",
                          "j"='NCI-H23',"c2"="CCLE"),
                     list("i"='REH', "c1"="GDSC",
                          "j"='REH',"c2"="CCLE"))
  transformed.list <- list(list("i"='KNS-81-FD', "c1"="GDSC",
                                "j"='KNS-81-FD',"c2"="CCLE"),
                           list("i"='COLO-320-HSR', "c1"="GDSC",
                                "j"='COLO-320-HSR',"c2"="CCLE"))
  validation.list <- list(list("i"='OCI-LY10', "c1"="CCLE",
                               "j"='KARPAS-422',"c2"="GDSC"))
  
  setwd(out.dir)
  discordant.x <- runAndPlotSummary(discordant.list, "discordant", 
                                    analysis.status='within', max.diff=1)
  drift.x <- runAndPlotSummary(drift.list, "drift", 
                               analysis.status='within', max.diff=1)
  transformed.x <- runAndPlotSummary(transformed.list, "transformed", 
                                     analysis.status='within', max.diff=1)
  validation.x <- runAndPlotSummary(validation.list, "validation", 
                                    analysis.status='within', max.diff=1)
}
