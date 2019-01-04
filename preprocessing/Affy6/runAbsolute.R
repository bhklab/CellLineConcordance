.libPaths(c(.libPaths(), "/mnt/work1/users/home2/quever/R/x86_64-unknown-linux-gnu-library/3.0"))
library(ABSOLUTE)

#### SUBFUNCTIONS 
# Name:       doAbsolute
# Purpose:    Sets the ABSOLUTE parameters and runs the RunAbsolute Command
# Arguments:  arg = absolute path to the HAPSEG file
#             d.name = cell line disease name
doAbsolute <- function(arg, d.name){
  #removes the .RData tag from the filenames if they exist
  p.name <- getSampleName(arg)
  #Create output Absolute directory path
  r.dir <- getAbsoluteDir(arg)
  
  #I/O Parameters
  hapseg.file <- arg  # Absolute filepath to Hapseg file
  out.dir <- r.dir  # Absolute filepath to an output directory (will be created)
  plot.name <- p.name   # name of the sample, included in output plots
  disease.name <- d.name  #  disease of the sample (used for display purpose)
  plat <- 'SNP_6.0'   # Supported values: ‘SNP_250K_STY’, 'SNP_6.0' and 'Illumina_WES'
  cn.type <- 'allelic'  # HAPSEG -> 'allelic || Segmentation -> 'total'
  verbose.val <- TRUE
    
  #ABSOLUTE Calculation Parameters
  max.sigma <- 0.02   # max value of excess sample level variance - Default:0.015
  sigma <- 0  # value of excess sample level variance used for mode search - Default:0
  
  min.pl <- 0.95  # min ploidy value to consider - Default:0.95
  max.pl <- 10  # max ploidy value to consider - Default:10
  
  max.seg <- 1500   # max number of allelic segments - Default: 1500
  max.neg <- 0  # max fraction of the genome modelled as less than zero - Default:0.005
  max.subclonal <- 0    # max fraction of genome that can be modelled as subclonal   - Default:0.05
  
  if (!file.exists(out.dir)) {   
    dir.create(out.dir, 
               recursive=TRUE)  
  }
  
  RunAbsolute(seg.dat.fn=hapseg.file,
              platform=plat,
              primary.disease=disease.name,
              sample.name=plot.name,
              results.dir=out.dir,
              copy_num_type=cn.type,
              
              max.sigma.h=max.sigma,
              sigma.p=sigma,
              min.ploidy=min.pl,
              max.ploidy=max.pl,
              max.as.seg.count=max.seg,
              max.neg.genome=max.neg,
              max.non.clonal=max.subclonal,
              
              verbose=verbose.val
              )
}

# Name:       parseHapList
# Purpose:    Parses the hap_list.txt containing absolute paths to HAPSEG output
#           .RData files
# Arguments: hap.path = Location of the HAPSEG file list
parseHapList <- function(hap.path){
  hap.df <- read.table(hap.path, header=TRUE, fill=TRUE)
  hap.vector <- hap.df$hap_files
  return(hap.vector)
}

# Name:       doReviewObject
# Purpose:    Combine all generated ABSOLUTE output files to run the function
#           CreateReviewObject for visualization and analysis purposes
# Arguments:  arg = vector containing absolute paths to the HAPSEG files
# Returns:    obj = String - descriptive name of collection of samples
doReviewObject <- function(arg){
  #removes the .RData tag from the filenames if they exist
  hap.name <- getSampleName(arg)
  abs.name <- sub("HAPSEG_", "", hap.name)
  abs.name <- sub("_segdat", "", abs.name)
  abs.name <- paste(abs.name, ".ABSOLUTE.RData", sep="")
  
  #Create output Absolute directory and Results path
  abs.dir <- getAbsoluteDir(arg)
  out.dir <- dirname(as.character(arg[1]))
  out.dir <- sub("\\/hapseg_output", "", out.dir)
  
  out.dir <- file.path(out.dir, "absolute", "summary")  #directory path to place the results
  #Compose a vector of all 
  absolute.f <- file.path(abs.dir, abs.name)  #vector of filenames pointing to ABSOLUTE output
  
  obj <- "ABS_summary"  #descriptive name for this collection of samples
  plot.ploidy <- TRUE   #FALSE to disable plotting of purity/ploidy modes
  plot.called <- FALSE  #TRUE to only plot purity/ploidy modes
  verbose.val <- TRUE   #verbose output
  cn.val <- 'allelic'   # HAPSEG -> 'allelic || Segmentation -> 'total'
  
  
  CreateReviewObject(obj.name=obj, 
                     absolute.files=absolute.f, 
                     indv.results.dir=out.dir, 
                     plot.modes=plot.ploidy,
                     copy_num_type='allelic',
                     verbose=verbose.val)
  return(obj)
}

# Name:       doExtractResults
# Purpose:    Produces the finalized output of ABSOLUTE CN caller
# Arguments:  arg = vector containing absolute paths to the HAPSEG files
#             obj = String - descriptive name of collection of samples
doExtractResults <- function(arg, obj){
  base.dir <- dirname(as.character(arg[1]))
  base.dir <- sub("\\/hapseg_output", "", base.dir)
  
  summary.dir <- file.path(base.dir, "absolute", "summary")  #directory path of summary
  output.dir <- file.path(base.dir, "absolute", "abs_extract") #directory path to place the results
  
  modes.path <- file.path(summary.dir, paste(obj, ".PP-modes.data.RData", sep=""))  #created from createReviewObject
  calls.path <- file.path(summary.dir, paste(obj, ".PP-calls_tab.txt", sep=""))   #name of the file to be uploaded
  analyst <- "BHK_lab"  #Analyst ID
  cn.type <- 'allelic'  # HAPSEG -> 'allelic || Segmentation -> 'total'
  
  ExtractReviewedResults(reviewed.pp.calls.fn=calls.path,
                         analyst.id=analyst,
                         modes.fn=modes.path,
                         out.dir.base=output.dir,
                         obj.name=obj,
                         copy_num_type=cn.type)
}

# Name:       getAbsoluteDir
# Purpose:    Converts absolute pathname to a file to absolute path to new ABSOLUTE dir
# Arguments:  arg = absolute path to the output HAPSEG file
# Returns:    r.dir = String  - absolute path to new directory to be created
getAbsoluteDir <- function(raw.dir){
  #Create output Absolute directory
  r.dir.raw <- dirname(as.character(raw.dir))
  r.dir.filt <- sub("\\/hapseg_output", "", r.dir.raw)
  r.dir <- file.path(r.dir.filt, "absolute", basename(as.character(raw.dir)))
  return(r.dir)
}

# Name:       getSampleName
# Purpose:    Converts absolute pathname to a file to absolute path to new ABSOLUTE dir
# Arguments:  arg = absolute path to the output HAPSEG file
# Returns:    p.name = String  - name of the sample
getSampleName <- function(raw.name){
  p.name <- sub("\\.RData$", "", raw.name)
  p.name <- basename(as.character(p.name))
  return(p.name)
}






#### MAIN 
args <- commandArgs(TRUE)
hap.file <- args[1]   # File containg all CEL files being analyzed
print(paste("HAPSEG output file list path: ", hap.file, sep=""))
hap.list <- parseHapList(hap.file)
print(paste("Hap List: ", hap.list, sep=""))

#results.dir <- '/mnt/work1/users/home2/quever/HK_lab/hapseg_output/'
#snp.fn <- '/mnt/work1/users/home2/quever/HK_lab/affy_output/quant-norm.pm-only.med-polish.expr.summary.txt'
#calls.fn <- '/mnt/work1/users/home2/quever/HK_lab/birdseed/birdseed-v1.calls.txt'
#clusters.fn <- '/mnt/work1/users/home2/quever/HK_lab/birdseed/birdseed-v1.snp-models.txt'
#cel.list <- parseCelList(cel.p)

for(hap in hap.list){
  print(paste("Analyzing ", basename(hap), sep=""))
  doAbsolute(hap, 'cooties')
}

obj <- doReviewObject(hap.list)

doExtractResults(hap.list, obj)
