.libPaths(c(.libPaths(), "/mnt/work1/users/pughlab/bin/r_lib/library/3.0"))
library(unixtools)
set.tempdir('/mnt/work1/users/pughlab/bin/r_lib/tmp')

library(foreach)
library(doParallel)

#### SUBFUNCTIONS 
# Name:       doHapseg
# Purpose:    Sets the HAPSEG parameters and runs the RunHapSeg Command
# Arguments:  arg = name of the .CEL file minus the .CEL extension
#             r.dir = output results.dir
#             snp.fn.loc = snp intensity file location
#             calls.fn.loc = birdseed calls file location
#             clusters.fn.loc = birdseed clusters file location
doHapseg <- function(arg, r.dir, snp.fn.loc, calls.fn.loc, clusters.fn.loc){
  .libPaths(c(.libPaths(), "/mnt/work1/users/pughlab/bin/r_lib/library/3.0"))
  library(HAPSEG)
  #Naming parameters
  plate <- 'HAPSEG'  # plate name
  array <- basename(arg)  # Name of the chip that was run
  results.dir.val <- file.path(r.dir, array)
  array <- sub("^X", "", array) # Removes the X in front of numbers for SU2C

  #Output from APT:
  snp.fn.val <- snp.fn.loc
  calls.fn.val <- calls.fn.loc
  clusters.fn.val <- clusters.fn.loc
  seg.fn.val <- NULL
  
  #HapSeg Parameters:
  genome <- 'hg19'  #Genome build - hg18 or hg19
  platform.val <- 'SNP_6.0' #Platform usage: SNP_250k_STY or SNP_6.0
  hapmap.pop <- 'CEPH' #CEPH, YOR, CH, JA (european, nigerian, chinese, japanese)
  # Phased BEAGLE files: ftp://ftp.broadinstitute.org/pub/genepattern/example_files/HAPSEG_1.1.1/phasedBGL.zip
  beagle.dir <- '/mnt/work1/users/pughlab/references/BEAGLE/phasedBGL'
  beagle.dir.genome <- file.path(beagle.dir, genome)
  
  #Segmentation Parameters:
  merge.threshold <- 0.0000000001 #Segmentation merge threshold - Default: 1e-10
  minimum.seg.size <- 5   # Minimum segment size - Default: 10
  outlier.prob <- 0.001  #outlier probability - Default:0.05
  
  #Default Set Parameters:
  # merge.small=TRUE    Merge small segments
  # merge.close=TRUE    Merge nearby segments
  # plot.segfit=TRUE    Plot segmentation export to jpg
  # drop.x=TRUE   Drop x-chrom from calculation
  # drop.y=TRUE   Drop y-chrom from calculation
  # verbose=TRUE  Be verbose
  # calibrate.data=TRUE   Calibration of input data using birdseed clusters file
  # imput.gt=TRUE   will impute genotypes via BEAGLE
  # normal=FALSE    Treat this sample as a normal
  # use.normal=TRUE   Use a matched normal if one is provided
  # snp.file.parser=AptSnpFileParser    Function used to parse the snp.file
  # clusters.file.parser=BirdseedClustersFileParser   Function used to parse the clusters.fn
  
  
  if (!file.exists(results.dir.val)) {   
    dir.create(results.dir.val, 
               recursive=TRUE)  
  }
  
  array <- gsub("^X", "", array)
  print(paste("Array name is: ", array, sep=""))
  sink(file=file.path(results.dir.val, paste(array, ".hapseg.out.txt", sep="")))
  RunHapSeg(plate.name=plate, 
            array.name=array,
            
            snp.fn=snp.fn.val, 
            calls.fn=calls.fn.val,
            clusters.fn=clusters.fn.val, 
            seg.fn=seg.fn.val, 
            calibrate.data=TRUE,
            genome.build=genome,
            platform=platform.val,  
            snp.file.parser=AptSnpFileParser, 
            clusters.file.parser=BirdseedClustersFileParser,
            
            use.pop=hapmap.pop,
            impute.gt=TRUE,
            phased.bgl.dir=beagle.dir.genome,
            normal=FALSE,
            use.normal=TRUE,
            
            min.seg.size=minimum.seg.size,
            seg.merge.thresh=merge.threshold,
            out.p=outlier.prob,
            merge.small=TRUE, 
            merge.close=TRUE, 
            plot.segfit=TRUE, 
            drop.x=TRUE,
            drop.y=TRUE,
            
            results.dir=results.dir.val, 
            
            adj.atten=FALSE,
            verbose=TRUE)
    sink()
}

# Name:       parseCelList
# Purpose:    Parses the .CEL list created by the perl script to create a 
#             vector containing all .CEL files
# Arguments: cel.path = Location of the .CEL file list
parseCelList <- function(cel.path){
  cel.df <- read.table(cel.path, header=TRUE, fill=TRUE, check.names=FALSE)
  #removes the .CEL tag from the filenames if they exist
  cel.vector <- sub("\\.CEL$", "", cel.df$cel_files, ignore.case=TRUE)
#  cel.vector <- sub("\\.CEL$", "", cel.df$cel_files)
  return(cel.vector)
}




#### MAIN 
args <- commandArgs(TRUE)
        
#Set up the clusters for parallel processing
num.of.cores <- 1

print(paste("Number of cores allocated: ", num.of.cores, sep=""))
#clust <- makeCluster(num.of.cores, outfile="")
#registerDoParallel(clust)

cel.p <- args[1]   # File containg all CEL files being analyzed
results.dir <- args[2]   # Output results directory
snp.fn <- args[3]    #snp intensity file location           
calls.fn <- args[4]    #calls.fn.loc = birdseed calls file location
clusters.fn <- args[5]   #clusters.fn.loc = birdseed clusters file location
print(paste("Cel Path: ", cel.p, sep=""))
cel.list <- parseCelList(cel.p)
print(paste("Cel List: ", cel.list, sep=""))

#results.dir <- '/mnt/work1/users/home2/quever/HK_lab/hapseg_output/'
#snp.fn <- '/mnt/work1/users/home2/quever/HK_lab/affy_output/quant-norm.pm-only.med-polish.expr.summary.txt'
#calls.fn <- '/mnt/work1/users/home2/quever/HK_lab/birdseed/birdseed-v1.calls.txt'
#clusters.fn <- '/mnt/work1/users/home2/quever/HK_lab/birdseed/birdseed-v1.snp-models.txt'
#cel.list <- parseCelList(cel.p)

for(cel in cel.list){
  print(paste("Analyzing ", cel, sep=""))
  x <- doHapseg(cel, results.dir, snp.fn, calls.fn, clusters.fn)
  print("Hapseg complete.  Saving results")
  #save(x, file=file.path(results.dir, paste(cel, ".R", sep="")))
}
