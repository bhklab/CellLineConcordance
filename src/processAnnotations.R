# Function: matchClAnno
# Purpose:  Given a threshold, it will go through each cell line and determine which cell lines match and which don't match
# Input:  x <- a row in cell.line.anno dataframe, that contains 'cgp.filename' and 'ccle.filename'
#         match_val <- A threshold value (default = 0.8 for snp, 0.6 for cnv)
#         matrix.perc.match <-  a matrix containing all possibel pair-wise concordance scores
# Returns:  List containing matches and mismatches and which cell lines are involved
matchClAnno <- function(x, match_val, matrix.perc.match, fn.headers, anno.mismatch=FALSE, verbose=FALSE){
  if(verbose) cat(paste(x['unique.cellid'], " - ", paste(fn.headers, collapse="/"), "\n", sep=""))
  match.values <- c()
  
  cl.fn <- as.vector(x[fn.headers])
  cl.fn <- cl.fn[!is.na(cl.fn)]
  if(length(grep(".cel", colnames(matrix.perc.match)[1], ignore.case=TRUE)) > 0) {
    if(length(grep(".cel", cl.fn, ignore.case=TRUE)) < length(cl.fn)){
      cl.fn <- as.vector(sapply(cl.fn, function(y) {
        if(length(grep(".cel", y, ignore.case=TRUE)) == 0)  y <- paste(y, "CEL", sep=".")
        return(y)
      }))
      
    }
  } else {
    cl.fn <- gsub(".cel", "", cl.fn, ignore.case = TRUE)
  }
  
  rm.fn <- vector(mode="character", length=0)

  # Removes cell lines in matrix.perc.match that are not annotated
  match.idx <- sapply(cl.fn, function(x) grep(paste("^", x, "((\\.hscr|\\.cel)|$)", sep=""), 
                                              colnames(matrix.perc.match), 
                                              ignore.case=TRUE, perl=TRUE))
  match.stat <- sapply(match.idx, length)
  if(any(match.stat == 0L)){
    rm.fn <- cl.fn[match.stat == 0L]
  } 
  
  # Sets a state no CEL files were found for the given cell line
  if(length(rm.fn) == length(cl.fn)){
    match.values <- -1
    cl.fn <- 'NA'
    
  # Finds all the matching or mismatching based on the annotated matching cl.fn list
  } else {
    if(match_val == "match"){
      cl.fn <- cl.fn[match.stat == 1L]
      match.values <- matrix.perc.match[cl.fn, cl.fn, drop=FALSE]
      match.values <- setAnnoNames(match.values)
    } else if(match_val == "mismatch"){
      cl.fn <- cl.fn[match.stat == 1L]
      mismatch.idx <- as.integer(match.idx[match.stat == 1L])
      cl.mismatch.fn <- colnames(matrix.perc.match)[-mismatch.idx]
      match.values <- matrix.perc.match[cl.fn, cl.mismatch.fn, drop=FALSE]
      na.idx <- which(apply(match.values,2, is.na))
      if(length(na.idx) > 0) match.values <- match.values[,-na.idx]
      if(anno.mismatch) match.values <- setAnnoNames(match.values)
    } 
  }
  
  # Catch-statement if the CEL file was not found in matrix.perc.match
  if(class(match.values) == "numeric"){
    #print(paste("No matching CEL files for: ", x['unique.cellid'], sep=""))
    match.values <- t(match.values)
  }
  return(as.matrix(match.values))
}


setAnnoNames <- function(match.values){
  getAnno <- function(cel.names){
    cel.idx <- lapply(cel.names, function(x) which(paste(x, ".CEL", sep="") == cell.line.anno, 
                                               arr.ind=TRUE))
    if(!class(cel.idx) %in% 'list') cel.idx <- list(cel.idx)
    ds.idx <- sapply(cel.idx, function(x) x[2])
    uid.idx <- sapply(cel.idx, function(x) x[1])
    fn.hd <- gsub(".filename(.+)?", "", colnames(cell.line.anno)[ds.idx], perl=TRUE)
    fn.hd <- paste(fn.hd, cell.line.anno[uid.idx, 'unique.cellid'], sep="_")
    fn.hd[which(is.na(ds.idx))] <- cel.names[which(is.na(ds.idx))]
    
    return(fn.hd)
  }
  
  
  colnames(match.values) <- getAnno(colnames(match.values))
  rownames(match.values) <- getAnno(rownames(match.values))
  return(match.values)
}


removeExtraCol <- function(x, row.name, match.stat){
  # Identify all non-(-1) matches
  x[is.na(x)] <- -1
  imp.idx <- which(apply(x, 2, function(y) length(y[unique(y) != -1]) >= 1))
  imp.values <- x[,imp.idx, drop=FALSE]

  
  if(dim(imp.values)[2] == 1){
    imp.headers <- colnames(x)[apply(x, 2, function(y) length(y[unique(y) != -1]) >= 1)]
    if(length(imp.headers) == 1){
      colnames(imp.values) <- imp.headers
      if(dim(imp.values)[1] == 1){
        rownames(imp.values) <- row.name
      }
    } else {
      imp.values <- t(imp.values)
      if(dim(imp.values)[1] == 1){
        rownames(imp.values) <- row.name
      }
    }
    
    if(dim(imp.values)[1] == 1){
      rownames(imp.values) <- row.name
    }
  } 
  
  if(match.stat == 'mismatch' & dim(imp.values)[2] > 1) imp.values <- setAnnoNames(imp.values)
  
  return(imp.values)
}