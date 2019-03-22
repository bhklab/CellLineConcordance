plotCytobandDiff <- function(diff.list, chr.df, max.diff=5){
  cn.diff.df <- diff.list[['cytoband']]
  ccl.data <- diff.list[['ccl.data']]
  cn.diff.df$cumStart <- as.numeric(as.character(cn.diff.df$cumStart))
  cn.diff.df$cumEnd <- as.numeric(as.character(cn.diff.df$cumEnd))
  diff.frac <- round((as.numeric(as.character(cn.diff.df$diff.val)) / max.diff), 4)
  diff.frac[diff.frac > 1] <- 1
  cn.diff.df$diff.frac <- diff.frac
  
  split.screen(matrix(c(0, 0.3, 0, 1,
                        0.3, 1, 0, 1), byrow = TRUE, ncol=4))
  # Concordances screen
  screen(1)
  par(mar=c(5.1, 0.25, 4.1, 0.25))
  barplot(do.call("cbind", ccl.data[c("geno", "a1", "a2", "tcn")]), xlim=c(0,1), 
          cex.axis = 0.5, horiz = TRUE)
  
  
  # CNV differences between both samples screen
  screen(2)
  par(mar=c(5.1, 0.25, 4.1, 0.25))
  plot(0, type='n', xlim=c(0, max(chr.df$cumEnd)), ylim=c(0,1), xlab='', ylab='', yaxt='n', xaxt='n')
  apply(cn.diff.df, 1, function(each.row){
    rect(xleft=each.row['cumStart'], ybottom = 0,
         xright = each.row['cumEnd'], ytop = 1, 
         col = alpha("red", each.row['diff.frac']), border = FALSE, lwd=0)
  })
  abline(v = chr.df$cumStart, lty=2, col="grey")
  axis(side = 1, at = chr.df$cumStart + with(chr.df, (cumEnd-cumStart)/2),
       labels = chr.df$chr, tick = FALSE, cex.axis=0.5)
  close.screen(all.screens=TRUE)
}

runAndPlotSummary <- function(list.x, file.grp.id='transformed', max.diff=4, ...){
  # Get Summary Info
  x.diff <- lapply(list.x, runListFun, ...)
  names(x.diff) <- sapply(list.x, getDsName)
  
  # Plot Summary information
  pdf(file.path("cnv_summary", paste0(file.grp.id, ".pdf")), height=2.5)
  lapply(x.diff, function(diff.list) plotCytobandDiff(diff.list, chr.df, max.diff = max.diff))
  dev.off()
  
  cat(paste("scp quever@mordor:", 
            file.path(getwd(), "cnv_summary", paste0(file.grp.id, ".pdf")), " .\n", sep=""))
  
  return(x.diff)
}

# Obtains the naming convention for ds.name
getDsName <- function(x){
  ds.name <- paste("all", paste(x$c1, x$i, sep="-"), 
                   paste(x$c2, x$j, sep="-"), sep=".")
  return(ds.name)
}

getGenoConc <- function(x, i, j, c1, c2){
  id1 <- cell.line.anno[which(cell.line.anno$unique.cellid == i), grep(paste0("^", c1, "\\.filename"), 
                                                                       colnames(cell.line.anno),
                                                                       ignore.case = TRUE)[1]]
  id2 <- cell.line.anno[which(cell.line.anno$unique.cellid == j), grep(paste0("^", c2, "\\.filename"), 
                                                                       colnames(cell.line.anno), 
                                                                       ignore.case = TRUE)[1]]
  geno.conc <- tryCatch({matrix.perc.match[id1, id2]}, error=function(e){NULL})
  names(x) <- c("a1", "a2", "tcn")
  append(x, list('geno'=geno.conc))
}

# Runs a list of cell lines throuhg the pipeline
runListFun <- function(x, cor.stat=TRUE, ...){
  i=x$i
  j=x$j
  c1=tolower(x$c1)
  c2=tolower(x$c2)
  ds.name <- paste("all", paste(c1, i, sep="-"), paste(c2, j, sep="-"), sep=".")
  
  adj.diff.data <- plotAllHscr(i=i, j=j, 
                               c1=c1, c2=c2, cor.stat=cor.stat, ...)
  if(cor.stat) adj.diff.data[['ccl.data']] <- getGenoConc(adj.diff.data[['ccl.data']], i, j, c1, c2)
  
  save(adj.diff.data,i, j, c1, c2, 
       file=file.path("adj_diff_data", paste(ds.name, ".adjDiff.rdata", sep="")))
  cytoband.diff <- processCytobandDiff(i=i, j=j, c1=c1, c2=c2, adj.diff.data=adj.diff.data)
  
  return(list("cytoband"=cytoband.diff,
              "ccl.data"=adj.diff.data[['ccl.data']]))
}

# Obtains and post-processes the cytoband IDs for each segment that was detected as different
processCytobandDiff <- function(i, j, c1, c2, adj.diff.data=adj.diff.data, ...){
  ds.name <- paste("all", paste(c1, i, sep="-"), paste(c2, j, sep="-"), sep=".")
  out.file <- file.path(getwd(), "cnv_cytoband", paste(ds.name, ".cytoband.pdf", sep=""))
  out.tsv <- file.path(getwd(), "cnv_cytoband", paste(ds.name, ".cytoband.tsv", sep=""))
  
  seg.diff <- getCytoband(ccl.diff=adj.diff.data[['ccl.diff']], ord.df, 
                          quant.thresh=0.9, min.thresh=0.2,
                          visualize = TRUE, 
                          chr.changepoints = adj.diff.data[['chr.changepoint']],
                          outfile=out.file)
  
  seg.diff <- collapseSegs(seg.diff, min.dist=10000)
  seg.diff <- filtSegs(seg.diff, min.seg.size=20003)
  seg.diff$seg.length <- round(seg.diff$seg.length / 1e6, 1)
  
  seg.diff <- cumSumCytoband(seg.diff)
  
  write.table(seg.diff, file = out.tsv, sep="\t", 
              quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  cat(paste("scp quever@mordor:", out.file, " .\n", sep=""))
  cat(paste("scp quever@mordor:", out.tsv, " .\n", sep=""))
  return(seg.diff)
}

# Takes the seg.diff or cytoband.diff datastructure and reports the sgment start and
# end position as a cumulative start/end for easier plotting
cumSumCytoband <- function(x){
  x <- apply(x, 2, function(col.x){
    gsub(" ", "", col.x)
  })
  
  # Obtains the cumulative start point for each chr and uses that to create
  # the cumulative start and end for each segment
  cum.x <- apply(x, 1, function(row.x) {
    cumStart.val <- chr.df[which(as.character(row.x['chr']) == chr.df$chr), 'cumStart']
    cumStart <- cumStart.val + as.integer(as.character(row.x['start.seg']))
    cumEnd <- cumStart.val + as.integer(as.character(row.x['end.seg']))
    return(c(row.x, "cumStart"=cumStart, "cumEnd"=cumEnd))
  })
  cum.x <- as.data.frame(t(cum.x))
  return(cum.x)
}

## Collapse consecutive segments that are a minimum distance apart
collapseSegs <- function(seg, min.dist=5000){
  # Identifies all contiguous segments that are within a set minimum distance
  collapse.idx <- sapply(c(2:nrow(seg)), function(x){
    if(seg[(x-1),]$chr == seg[x,]$chr) {
      pre.seg.val <- abs(seg[(x-1),]$end.seg - seg[x,]$start.seg)
    } else {
      pre.seg.val <- 1e9
    }
    pre.seg.val <= min.dist
  })
  
  # Identifies changepoints
  collapse.idx <- which(collapse.idx)
  collapse.diff.idx <- which(diff(collapse.idx) != 1)
  
  collapse.list <- list()
  start.idx <- 1
  # Constructs a "collapsed" row by merging the first and last segment in a contiguous segment
  for(cidx in collapse.diff.idx){
    end.idx <- cidx
    start.row <- seg[collapse.idx[start.idx],]
    end.row <- seg[collapse.idx[end.idx]+1,]
    
    t.seg.length <- sum(seg[c(collapse.idx[start.idx]:collapse.idx[end.idx]),]$seg.length)
    x.df <- data.frame("chr"= end.row$chr,
                       "start.seg"= start.row$start.seg,
                       "end.seg" = end.row$end.seg,
                       "diff.val" = with(seg[c(collapse.idx[start.idx]:collapse.idx[end.idx]),], 
                                         weighted.mean(diff.val, seg.length/t.seg.length)),
                       "seg.length" = end.row$end.seg - start.row$start.seg,
                       "seg.cytoband.val"=NA)
    start.cytoband <- gsub("-.*", "", start.row$seg.cytoband.val)
    end.cytoband <- gsub("chr.*:(.*-)?", "", end.row$seg.cytoband.val)
    if(length(grep(end.cytoband, start.cytoband)) > 0){
      x.df$seg.cytoband.val <- start.cytoband
    } else {
      x.df$seg.cytoband.val <- paste(start.cytoband, end.cytoband, sep="-")
    }
    collapse.list [[as.character(cidx)]] <- x.df
    
    start.idx <- end.idx + 1
  }
  
  # Reconstructs the original segment df to remove the hyper-segmented segments and re-order
  seg2 <- rbind(seg[-collapse.idx,],
                do.call("rbind", collapse.list))
  seg <- seg2[with(seg2, order(chr, start.seg)),]
  seg
}



## Filters out segments that are smaller than a set size
filtSegs <- function(seg, min.seg.size = 20003){
  small.window.idx <- which(seg$seg.length <= min.seg.size)
  
  # Check which windows are non-contiguous (i.e. small windows that are segmented due to changing adj.diff values)
  small.window.filt <- sapply(small.window.idx, function(x){
    pre.seg.val <- post.seg.val <- 2
    if(x > 1){
      if(seg[(x-1),]$chr == seg[x,]$chr) {
        pre.seg.val <- abs(seg[(x-1),]$end.seg - seg[x,]$start.seg)
      }
    }
    if(x < nrow(seg)){
      if(seg[(x+1),]$chr == seg[x,]$chr) {
        post.seg.val <- abs(seg[(x+1),]$start.seg - seg[x,]$end.seg)
      }
    }
    
    if(pre.seg.val == 1 | post.seg.val == 1) FALSE else TRUE
  })
  seg[-small.window.idx[small.window.filt],]
}


# Uses rle() to obtain consecutive values in a dataframe and index them
getRleIdx <- function(x, col.id, na.val=-100){
  
  if(length(col.id) > 1) {
    uniq.id <- apply(x, 1, function(y) paste(y[col.id], collapse="-"))
    x$uniq <- uniq.id
    col.id <- 'uniq'
  }
  reformat.na.x <- as.character(x[,col.id])
  reformat.na.x[which(is.na(reformat.na.x))] <- na.val
  rle.x <- rle(reformat.na.x)
  
  
  #Get the array index for the start-to-end of each unique value/changepoint
  rle.x$start.idx <- c(1, (cumsum(rle.x$lengths) + 1)[-length(rle.x$lengths)])
  rle.x$end.idx <- rle.x$start.idx + (rle.x$lengths - 1)
  rle.x$values[which(rle.x$values == na.val)] <- NA
  rle.x$na.stat <- !is.na(rle.x$values)
  return(rle.x)
}

scaledSdNorm <- function(dat, targ.sd){
  apply(dat, 2, function(x){
    ((x - mean(x))/(1/targ.sd))/(sd(x))
  })
}

meanCenter <- function(dat, targ.mean){
  apply(dat, 2, function(x) (x - mean(x)) + targ.mean)
}

getCytoband <- function(ccl.diff, ord.df, quant.thresh=0.9, min.thresh=0.2,
                        visualize=FALSE, chr.changepoints=NA, outfile=NA){
  ## Function: Takes the adj.diff (adjusted differences between CN-profiles)
  ##    and annotates with cytobands of differences that exceed a set threshold
  
  # Segments the CN differences into continuous segments
  ccl.diff.rle <- getRleIdx(as.data.frame(ccl.diff), 1)
  start.row <- ord.df[ccl.diff.rle$start.idx, ]
  end.row <- ord.df[ccl.diff.rle$end.idx, ]
  
  # Populates the segment cytoband dataframe
  seg.cytoband <- cbind(start.row[,c("Chr", "Start.bin")], 
                        end.row[,"End.bin"], 
                        as.numeric(ccl.diff.rle$values))
  rownames(seg.cytoband) <- ccl.diff.rle$start.idx
  colnames(seg.cytoband) <- c("chr", "start.seg", "end.seg", "diff.val")
  
  # Uses the CollapsABEL::cytoband() function to annotate segments
  seg.cytoband.val <- apply(seg.cytoband, 1, function(x){
    chr.val <- paste0("chr",  as.integer(x['chr']))
    start.idx <- cytoband(chr = chr.val,pos = as.integer(x['start.seg']),  ref = 'hg19')
    end.idx <- cytoband(chr = chr.val, pos = as.integer(x['end.seg']),  ref = 'hg19')
    if(start.idx != end.idx) {
      cytoband.idx <- paste(chr.val, paste(c(start.idx, end.idx), collapse="-"), sep=":")
    } else {
      cytoband.idx <- paste(chr.val, start.idx, sep=":")
    }
    return(cytoband.idx)
  })
  
  # Isolates segments to report that pass a cutoff threshold
  seg.cytoband <- cbind(seg.cytoband, 
                        "seg.length"=with(seg.cytoband, end.seg - start.seg), 
                        seg.cytoband.val)
  # diff.idx <- which(as.numeric(seg.cytoband$diff.val) > quantile(ccl.diff, quant.thresh) |
  #                     as.numeric(seg.cytoband$diff.val) > min.thresh)
  diff.idx <- which(as.numeric(seg.cytoband$diff.val) > min.thresh)
  seg.cytoband <- seg.cytoband[diff.idx,]
  
  #Visualize the barplot and collapsed RLE indexed segments to ensure they are the same
  #thresh.val <- max(min.thresh, quantile(ccl.diff, quant.thresh))
  thresh.val <- min.thresh
  if (visualize) visualizeCytobands(ccl.diff, outfile, chr.changepoints,
                                    ccl.diff.rle, seg.cytoband, thresh.val)
  return(seg.cytoband)
}


#Visualize Cytobands
visualizeCytobands <- function(ccl.diff, outfile, chr.changepoints,
                               ccl.diff.rle, seg.cytoband, thresh.val){
  ## Plotting the barplot of adj.diff as well as the rle lineplots to ensure they are the same
  pdf(outfile)
  split.screen(c(2,1))
  screen(1)
  par(mar=c(0.5, 4.1, 2.1, 2.1))
  barplot.idx <- barplot(ccl.diff, las=2)
  abline(v = barplot.idx[chr.changepoints], col="grey", lty=4)
  
  screen(2)
  par(mar=c(2.1, 4.1, 0.5, 2.1))
  plot(0, type='n', ylim=c(0, 0.6), xlim=c(0, length(ccl.diff)), xaxt='n', las=2)
  lines(x=as.integer(t(cbind(ccl.diff.rle$start.idx, ccl.diff.rle$end.idx))),
        y=as.numeric(t(cbind(ccl.diff.rle$values, ccl.diff.rle$values))),
        col="red")
  abline(v = chr.changepoints, col="grey", lty=4)
  abline(h = thresh.val, col="grey", lty=4)
  axis(side = 1, at = chr.changepoints, labels = seq(1, 22), cex.axis=0.5, adj=1)
  text(x = rownames(seg.cytoband), y=0.4, 
       labels=seg.cytoband$seg.cytoband.val, adj = 0, cex=0.4, srt=90)
  close.screen(all.screens=TRUE)
  dev.off()
}

#compareCn function for squared-difference between cn-bins
compareCn <- function(cn1, cn2, median.norm=0, mean.norm=0, quant.norm=0, ret='sd'){
  if(median.norm==1){
    print("Median centering the segments")
    avg.median <- mean(c(median(cn1), median(cn2)))
    
    cn1 <- ((cn1 - median(cn1)) + avg.median)
    cn2 <- ((cn2 - median(cn2)) + avg.median)
    
  } else if (mean.norm == 1){
    print("Mean-centering, variance normalizing the segments")
    avg.mean <- mean(c(mean(cn1), mean(cn2)))
    avg.sd <- mean(c(sd(cn1), sd(cn2)))
    
    mat.cn <- matrix(c(cn1, cn2), byrow=FALSE, ncol=2)
    mat.cn.norm <- scaledSdNorm(mat.cn, avg.sd)
    mat.cn.norm.meanCenter <- meanCenter(mat.cn.norm, avg.mean)
    
    cn1 <- mat.cn.norm.meanCenter[,1]
    cn2 <- mat.cn.norm.meanCenter[,2]
  } else if(quant.norm == 1){
    print("Quantile normalizing the segments")
    require(preprocessCore)
    
    cn.mat <- matrix(c(cn1, cn2), ncol=2, byrow = FALSE)
    colnames(cn.mat) <- c("cn1", "cn2")
    cn.qmat <- normalize.quantiles(cn.mat, copy=TRUE)
    colnames(cn.qmat) <- c("cn1", "cn2")
    
    cn1 <- cn.qmat[,'cn1']
    cn2 <- cn.qmat[,'cn2']
  }
  
  #SSD Value (in -log space, normalized to worst case scenario)
  if(ret %in% 'sd'){
    cn.squared.diff <- (cn1 - cn2)^2
  } else {
    list("cn1" = cn1, "cn2" = cn2)
  }
}

getCclUid <- function(uid){
  uid.idx <- grep(paste("^", uid, "$", sep=""), cell.line.anno$unique.cellid, ignore.case=TRUE)
  if(length(uid.idx) == 0) uid.idx <- which(cell.line.anno$unique.cellid %in% uid)
  cl.uid <- cell.line.anno[uid.idx,]$unique.cellid
  names(cl.uid) <- uid.idx
  return(cl.uid)
}

getCclFile <- function(uid, ds){
  ccl.file <- cell.line.anno[as.integer(names(uid)), fn.headers[tolower(ds)]]
  ccl.file <- gsub(".cel", "", ccl.file, ignore.case=TRUE)
  return(ccl.file)
}

printStatus <- function(hscr, cl1, ds1, file1, cl2, ds2, file2){
  print(paste("Processing ", hscr, sep=""))
  
  print(paste("Cell Line UID: ", cl1, sep=""))
  print(paste("      dataset: ", ds1, sep=""))
  print(paste("         file: ", file1, sep=""))
  
  print(paste("Cell Line UID: ", cl2, sep=""))
  print(paste("      dataset: ", ds2, sep=""))
  print(paste("         file: ", file2, sep=""))
  print("#-------------")
}


loadStdev <- function(hscr){
  load(file.path(ccl.cnv.datadir, stdev.dir, paste(hscr, 'stdevMat.Rdata', sep=".")))
  names(stdev.list) <- toupper(names(stdev.list))
  return(stdev.list)
}

# Load in the raw CN data for each CEL file from each dataset
loadHscrMat <- function(path.to.dir, hscr, profile.id="mem.cnProfile", annofile=NULL){
  cwd <- getwd()
  setwd(path.to.dir)
  path.to.file <- file.path(path.to.dir, paste(hscr, profile.id, 'desc', sep="."))
  shared.desc <- dget(path.to.file)
  hscr.mat <- attach.big.matrix(shared.desc)
  
  options(bigmemory.allow.dimnames=TRUE)
  if(!is.null(annofile)){
    colnames(hscr.mat) <- mapIds(NULL, annofile, hscr.mat, NULL, 'colnames')
  } else {
    colnames(hscr.mat) <- gsub(".hscr.*Rdata", "", colnames(hscr.mat))
  }
  
  setwd(cwd)
  if(hscr == 'hscra1') {
    load(file.path(path.to.dir, "bin_order.RData"))
    return(list('hscr'=hscr.mat, 'ord.df'=ord.df))
  } else {
    return(hscr.mat)
  }
}

#' darken
#' @description https://gist.github.com/Jfortin1/72ef064469d1703c6b30
#' @param color Color value (e.g. "#fdc086")
#' @param factor Numeric value to scale the colour by [default=1.4]
#'
#' @return Color value that is darkened 
#' @export
#'
#' @examples
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- tryCatch({ 
    rgb(t(col), maxColorValue=255)
  }, error=function(e){
    col.idx <- which(col < 0)
    col[col.idx] <- 0
    col[c(1:3)[-col.idx][1]] <- 50
    rgb(t(col), maxColorValue=255)
  })
  col
}

#' lighten
#'
#' @param color Color value (e.g. "#fdc086")
#' @param factor Numeric value to scale the colour by [default=1.4]
#'
#' @return Color value that is lightened 
#' @export
#'
#' @examples
lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- tryCatch({ 
    rgb(t(col), maxColorValue=255)
  }, error=function(e){
    col.idx <- which(col > 255)
    col[col.idx] <- 255
    col[c(1:3)[-col.idx][1]] <- 50
    rgb(t(col), maxColorValue=255)
  })
  col
}


#' writeBigMatrix
#'
#' @param object A matrix to write into bigmatrix
#' @param uniq.id Outfile ID for current working directory: [uniq.id].[profile.id].bin
#' @param profile.id Outfile ID for current working directory: [uniq.id].[profile.id].bin
#'
#' @return NA
#' @export
#'
#' @examples
writeBigMatrix <- function(object, uniq.id, profile.id){
  write.table(object, file=paste(uniq.id, "allFiles.tsv", sep="."), 
              sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  obj.mat <- read.big.matrix(
    paste(uniq.id, "allFiles.tsv", sep="."), 
    sep="\t", has.row.names=TRUE, type ="double", header = TRUE, 
    backingfile = paste(uniq.id, profile.id, "bin", sep="."), 
    descriptorfile = paste(uniq.id, profile.id, "desc", sep=".")) 
}


#' mapIdsToFiles
#'
#' @param i Cell line ID 1 (e.g. KNS-81-FD)
#' @param j Cell line ID 2 (e.g. MCF7)
#' @param mm.list A list containing all cell line pairs in the reference dataset and whether they have the matching annotation, as well as the CNV concordance value
#'
#' @return
#' @export
#'
#' @examples
mapIdsToFiles <- function(i, j=NA, c1, c2, hscr, mm.list, verbose=FALSE){
  if(is.na(j)) j <- i
  print(paste("Evaluating ", i, "-", j, sep=""))
  mm.df <- mm.list[[grep(paste0("^", i, "$"), names(mm.list))]]
  if(verbose) print(mm.df)
  
  if(1==0){ # Deprecated
    # Samples with (matching annotation + below threshold)
    # Samples with (non-matching annotation + above threshold)
    match.low.idx <- which(mm.df$anno %in% 'Matching annotation' & 
                             mm.df$value < threshold)
    mismatch.hi.idx <- which(mm.df$anno %in% 'Different annotation' & 
                               mm.df$value >= threshold)
    all.idx <- c(match.low.idx, mismatch.hi.idx)
  }
  
  match.idx <- which(apply(mm.df, 1, function(x){
    all(c(any(grepl(paste0("^", c1, "_", i), x)),
          any(grepl(paste0("^", c2, "_", j), x))))
  }))
  if(length(match.idx) > 0){   
    each.idx <-  as.integer(match.idx)[1]
  } else {    #If no known match is found
    print("No known match is found, checking to make sure the files exist")
    each.idx <- nrow(mm.df) + 1
    mm.df <- rbind(mm.df, data.frame("X1" = paste(c1, i, sep="_"),
                                     "X2" = paste(c2, j, sep="_"),
                                     "value"=1,
                                     "anno"=NA,
                                     "fingerprint"="CNV"))
  }
  print(mm.df[each.idx,])
  # X1             X2     value                anno fingerprint
  # 342 gdsc_KNS-81-FD ccle_KNS-81-FD 0.9387356 Matching annotation         CNV
  
  parseIds <- function(each.idx, mm.df, x.col='X1'){
    cl.spl <- unlist(regmatches(as.character(mm.df[each.idx, x.col]),
                                regexpr("_", as.character(mm.df[each.idx, x.col])), 
                                invert = TRUE))
    cl.x <- getCclUid(cl.spl[2])
    if(length(cl.x) == 0) cl.x <- cl.spl[2]
    ds.x <- toupper(cl.spl[1])
    file.x <- getCclFile(cl.x, ds.x)
    ds.x <- gsub("[1-3]$", "", ds.x)
    
    list("ds"=ds.x, "cl"=cl.x, "file"=file.x)
  }
  ids1 <- parseIds(each.idx, mm.df, 'X1')
  ids2 <- parseIds(each.idx, mm.df, 'X2')
  ds.name <- paste(hscr, 
                   paste(ids1[['ds']], ids1[['cl']], sep="-"), 
                   paste(ids2[['ds']], ids2[['cl']], sep="-"), sep=".")
  
  printStatus(hscr, 
              ids1[['cl']], ids1[['ds']], ids1[['file']],
              ids2[['cl']], ids2[['ds']], ids2[['file']])
  return(list("ids1"=ids1, "ids2"=ids2, "dsname"=ds.name))
}


#' getMmDf
#'
#' @param i Cell line ID 1 (e.g. KNS-81-FD)
#' @param j Cell line ID 2 (e.g. MCF7)
#' @param c1 Dataset ID 1 (e.g. gdsc)
#' @param c2 Dataset ID 2 (e.g. ccle)
#' @param gen.report Boolean:  Generates a TSV of the adjusted difference, CN, and weights
#' @param hscr Character: For labelling purpose of Alleles (e.g. hscra1)
#' @param hscr.mat1 Matrix: CN-bins x Samples with sample IDs on column names
#' @param stdev.list1 List containing each individual dataset and StdDev values stored in a similar matrix
#' @param hscr.mat2 Different matrix related to j or c2 if not part of the original CCLE/GDSC/CGP/Pfizer [Default=NULL]
#' @param stdev.list2 Different list related to j or c2 if not part of the original CCLE/GDSC/CGP/Pfizer [Default=NULL]
#' @param output.pdf Boolean: Generates an external visualization or not
#' @param get.cor Boolean: Returns a correlation value
#' @param plot.type Character: Either 'png' or 'pdf'
#' @param analysis.status Either 'within' if c1 and c2 are in hscr.mat1, or 'between' if hscr.mat2 is not NULL
#'
#' @return
#' @export
#'
#' @examples
getMmDf <- function(i, j=NA, c1, c2, gen.report=TRUE, hscr=hscr,
                    hscr.mat1=hscr.mat.a1, stdev.list1=stdev.list.a1,
                    hscr.mat2=NULL, stdev.list2=NULL, 
                    analysis.status='within',
                    output.pdf=FALSE, get.cor=FALSE, plot.type='png'){
  ## Obtain proper file IDs and reference matrices
  validateIJ <- function(x){
    if(grepl("[^a-zA-Z0-9]", x)) x <- gsub("[^a-zA-Z0-9]", ".", x)
    x
  }
  i <- validateIJ(i)
  j <- validateIJ(j)
  
  all.ids <- mapIdsToFiles(i, j, c1, c2, hscr, mm.list, verbose=FALSE)
  ds1 <- all.ids[['ids1']][['ds']];     
  ds2 <- all.ids[['ids2']][['ds']];
  cl1 <- all.ids[['ids1']][['cl']];     
  cl2 <- all.ids[['ids2']][['cl']]
  ds.name <- all.ids[['dsname']]
  
  if(analysis.status == "within"){
    file1 <- all.ids[['ids1']][['file']]; 
    file2 <- all.ids[['ids2']][['file']];
    ds.name <- all.ids[['dsname']]
    
    hscr.mat2 <- hscr.mat1
    stdev.list2 <- stdev.list1
  } else if(analysis.status == 'between' & !is.null(hscr.mat2)){
    print(paste(cl1, ds1))
    print(paste(cl2, ds2))
    file1 <- all.ids[['ids1']][['file']];
    file2 <- colnames(hscr.mat2)[grep(paste0("^", j, "$"), 
                                      colnames(hscr.mat2))]
    
    stdev.list2 <- list()
    stdev.list2[[ds2]] <- matrix(rep(1, nrow(hscr.mat2)), ncol=1, dimnames = list(NULL, file2)) # Set no weights
    
  }
  
  ## Validation to ensure files 
  if(!all(sapply(list(file1, file2), function(x) length(x) > 0))) stop("Files were not found")
  if(!all(sapply(list(grep(paste("^", file1, "$", sep=""), colnames(hscr.mat1)), 
                      grep(paste("^", file2, "$", sep=""), colnames(hscr.mat2))),
                 function(x) length(x) > 0))) warning("Files were not found in column headers")
  print("Validated")
  
  
  ## CN and Stdev values
  rep1 <- hscr.mat1[,file1]
  rep2 <- hscr.mat2[,file2]
  rep.l <- compareCn(rep1, rep2, mean.norm=1, ret='plot')
  
  file1.stddev <- stdev.list1[[ds1]][,file1]
  file2.stddev <- stdev.list2[[ds2]][,file2]
  weight.f <- (file2.stddev * file1.stddev)   # StdDev weights
  adj.diff <- ((rep.l[['cn1']] - rep.l[['cn2']])^2 * weight.f)
  
  
  ###############################
  #     Visualization
  my_palette <- colorRampPalette(c("black", "white"))(n = 1000)
  weight_palette <- colorRampPalette(c("red", "red", "black"))(n = 1000)
  
  #close.screen(all.screens=TRUE)
  ds.name <- gsub("\\/", "-", ds.name)
  if(plot.type=='png') png(file.path("cnv_plots", paste(ds.name, ".png", sep="")), 
                           width=2000, height=1750, res=80)
  if(plot.type=='pdf') pdf(file.path("cnv_plots", paste(ds.name, ".pdf", sep="")), 
                           width=23, height=26)
  # pdf(paste(ds.name, ".pdf", sep=""))
  split.screen(matrix(c(0,1,0,0.1,
                        0,1,0.1,0.15,
                        0,1,0.15,0.2,
                        0,1,0.2,0.55,
                        0,1,0.55,0.9,
                        0,1,0.9,1), ncol=4, byrow=TRUE))
  chr.changepoints <- c(1,1+which(diff(ord.df$Chr)!=0))
  
  # Screen for the gradient legend indicating weight-factors:
  screen(1)
  par(mar=c(0, 4.1, 0, 2.1))
  plot(0, type="n", axes = FALSE, ylab="", xlab="")
  # Screen for the colour plot on weights 1
  screen(2)
  par(mar=c(0, 4.1, 0, 2.1))
  file1.col <- ceiling(file1.stddev * 1000)
  file1.col[which(file1.col < 0)] <- 0
  barplot(rep(1, length(file1.stddev)),  
          ylim=c(0,1), yaxt='n', ylab=paste(ds1, cl1, sep="-"),
          col=my_palette[file1.col], border=NA)
  
  # Screen for the colour plot on weights 2
  screen(3)
  par(mar=c(0, 4.1, 0, 2.1))
  file2.col <- ceiling(file2.stddev * 1000)
  file2.col[which(file2.col < 0)] <- 0
  barplot(rep(1, length(file2.stddev)),  
          ylim=c(0,1), yaxt='n', ylab=paste(ds2, cl2, sep="-"),
          col=my_palette[file2.col], border=NA)
  
  
  # Screen for the raw Weight-adjusted difference barplots
  screen(4)
  par(mar=c(0, 4.1, 0, 2.1))
  weight.col <- ceiling(weight.f * 1000)
  
  barplot.idx <- barplot(adj.diff, ylim=c(0,2), col=weight_palette[weight.col],
                         border=NA)
  legend("topright", paste("SSD: ", sum(adj.diff), sep=""), bty="n")
  abline(v = 0)
  abline(v = barplot.idx[chr.changepoints], col="grey", lty=4)
  
  # Screen for the raw Allele specific plots
  screen(5)
  par(mar=c(1, 4.1, 0, 2.1))
  cn1.ds.col <- ds.col[toupper(ds1)]
  cn2.ds.col <- ds.col[toupper(ds2)]
  if(ds1 == ds2) cn2.ds.col <- lighten(cn2.ds.col)
  plot(0, type='n', xlim=c(0, length(rep.l[['cn1']])), 
       col=cn1.ds.col, xaxt='n', ylab="hscr", ylim=c(0,10), pch=15)
  points(rep.l[['cn1']], col=cn1.ds.col, pch=15)
  points(rep.l[['cn2']], col=cn2.ds.col, pch=15)
  legend("topright", c(paste(ds1, file1, sep="-"),
                       paste(ds2, file2, sep="-")), 
         lwd=c(2,2), lty=c(1,1), col=c(cn1.ds.col, cn2.ds.col))
  abline(v = chr.changepoints, col="grey", lty=4)
  text(x=(chr.changepoints + 500), y=9, labels=c(1:22), pos=4, cex=1)
  
  close.screen(all.screens=TRUE)
  dev.off()
  cat(paste("scp quever@mordor:", 
            file.path(getwd(), "cnv_plots", paste(ds.name, 
                                                  if(plot.type=='png') ".png" else ".pdf", sep="")), 
            " .\n", sep=""))
  
  if(gen.report){
    dir.create(file.path("cnv_tsv"), showWarnings = FALSE)
    write.table(x = cbind(ord.df, 
                          data.frame("adjusted_diff"=adj.diff,
                                     "cn1"=rep.l[['cn1']],
                                     "cn2"=rep.l[['cn2']],
                                     "wf"=weight.f)),
                file=file.path("cnv_tsv", paste(ds.name, ".tsv", sep="")), 
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    cat(paste("scp quever@mordor:", 
              file.path(getwd(), "cnv_tsv", paste(ds.name, ".tsv", sep="")), " .\n", sep=""))
  }
  r.grob <- c()
  if(output.pdf){
    require(gridExtra)
    require(png)
    print("Rasterizing png file...")
    r.grob <- rasterGrob(readPNG(file.path("cnv_plots", paste(ds.name, ".png", sep="")), native = FALSE),
                         interpolate = FALSE)
  }
  
  cn.wtd.cor <- NULL
  if (get.cor) cn.wtd.cor <- wtd.cor(rep1, rep2, file1.stddev * file2.stddev)[1]
  return(list("grob"=r.grob,
              "ccl.diff"=adj.diff,
              "chr.changepoint"=chr.changepoints,
              "wtd.cor"=cn.wtd.cor))
}

#' plotAllHscr
#'
#' @param i Cell line ID 1 (e.g. KNS-81-FD)
#' @param j Cell line ID 2 (e.g. MCF7)
#' @param c1 Dataset ID 1 (e.g. gdsc)
#' @param c2 Dataset ID 2 (e.g. ccle)
#' @param cor.stat Boolean: Returns a correlation value
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plotAllHscr <- function(i, j, c1, c2, cor.stat=FALSE, analysis.status='within', ...){
  # Error in getMmDf(i = i, j = j, c1 = c1, c2 = c2, hscr.mat = hscr.mat.a1,  : 
  #                    argument 5 matches multiple formal arguments
  print(analysis.status)
  h1 <- getMmDf(i = i, j = j, c1 = c1, c2 = c2, 
                hscr.mat1=hscr.mat.a1, stdev.list1=stdev.list.a1, hscr='hscra1', 
                output.pdf = TRUE, get.cor=cor.stat,
                hscr.mat2= if(analysis.status=='between') x.hscr.mat.a1 else NULL,
                analysis.status=analysis.status,
                ...)
  h2 <- getMmDf(i = i, j = j, c1 = c1, c2 = c2, 
                hscr.mat1=hscr.mat.a2, stdev.list1=stdev.list.a2, hscr='hscra2', 
                output.pdf = TRUE, get.cor=cor.stat,
                hscr.mat2= if(analysis.status=='between') x.hscr.mat.a2 else NULL,
                analysis.status=analysis.status,
                ...)
  h12 <- getMmDf(i = i, j = j, c1 = c1, c2 = c2, 
                 hscr.mat1=hscr.mat.a12, stdev.list1=stdev.list.a12, hscr='hscra12', 
                 output.pdf = TRUE, get.cor=cor.stat,
                 hscr.mat2= if(analysis.status=='between') x.hscr.mat.a12 else NULL,
                 analysis.status=analysis.status,
                 ...)
  ds.name <- paste("all", paste(c1, i, sep="-"), paste(c2, j, sep="-"), sep=".")
  
  require(ggplot2)
  pdf(file.path("cnv_plots", paste(ds.name, ".pdf", sep="")))
  lapply(list(h1[['grob']], h2[['grob']], h12[['grob']]), function(x){
    print(qplot(1:10, 1:10, geom="blank") +
            annotation_custom(x, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"),
                  axis.title.x=element_blank(), axis.title.y=element_blank(),
                  axis.text.x=element_blank(), axis.text.y=element_blank(),
                  axis.ticks.x=element_blank(), axis.ticks.y=element_blank()))
  })
  #do.call(grid.arrange, c(list(h1, h2, h12), nrow = 1, ncol=1))
  dev.off()
  cat(paste("scp quever@mordor:", 
            file.path(getwd(), "cnv_plots", paste(ds.name, ".pdf", sep="")), " .\n", sep=""))
  
  return(list("ccl.diff"=h12[['ccl.diff']],
              "chr.changepoint"=h12[['chr.changepoint']],
              "ccl.data"=list("h1"=h1[['wtd.cor']],
                              "h2"=h2[['wtd.cor']],
                              "h12"=h12[['wtd.cor']])))
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
      
      if(nrow(match.idx) == 0){
        col.idx <- grep(x, cell.line.anno)
        if(length(col.idx) > 0){
          row.idx <- grep(x, cell.line.anno[,col.idx])
          match.idx <- matrix(c(row.idx, col.idx), nrow=1)
        } else {
          print(paste0("Can not find ", x))
        }
      }
      
      if(nrow(match.idx)>1) match.idx <- match.idx[1, ,drop=FALSE]
      cell.id <- cell.line.anno[match.idx[,1],'unique.cellid']
      paste0(gsub(".filename.*", "", colnames(cell.line.anno)[match.idx[,2]]), 
             "_", cell.id)
    })
  }
  map.ids
}

