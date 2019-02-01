library(plyr)
library(magrittr)


setwd("/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2017_GNE/snp_fingerprint/output")
snp.conc.files <- list.files(pattern = "Rdata")


initializeConcordanceMatrix <- function(type="match"){
  if(type=='match'){
    datasets <- c("ccle", "gdsc", "cgp", "pfizer")
    geno.mat <- as.data.frame(matrix(ncol=length(datasets), nrow=0))
    colnames(geno.mat) <- datasets
  } else {
    geno.mat <- as.data.frame(matrix(ncol=2, nrow=0))
    colnames(geno.mat) <- c("stat", "ids")
  }
  geno.mat
}

interactive <- FALSE

match.anno.df <- initializeConcordanceMatrix()
nonmatch.anno.df <- initializeConcordanceMatrix(type="nonmatch")
#for(each.file in snp.conc.files){
{  
  print(each.file)
  load(each.file)
  
  if(!exists("top.conc")) {
    top.conc <- top.corr
    names(top.conc) <- colnames(corr.mat)
  }
  all.geno.mats <- lapply(top.conc, function(each.cl){
    print(colnames(each.cl)[1])
    geno.mat <- initializeConcordanceMatrix()
    matchCl <- function(df){
      ref.id <- gsub("_.+", "", rownames(df) )
      df <- matrix(df[,1], nrow=1, 
                   dimnames = list(colnames(df)[1],
                                   ref.id))
      df
    }
    
    checkMatches <- function(each.cl){
      rows <- gsub("^.+?_", "", rownames(each.cl)) %>%
        gsub("[^a-zA-Z0-9]", "", .) %>%
        gsub("[oO]", "0", .)
      cols <- gsub("[^a-zA-Z0-9]", "", colnames(each.cl)[1]) %>%
        gsub("[oO]", "0", .)
      
      m.idx <- sapply(rows, function(x) grepl(x, cols, ignore.case = TRUE))
      m.idx2 <- as.vector(sapply(cols, function(x) grepl(x, rows, ignore.case = TRUE)))
      idx <- apply(rbind(m.idx, m.idx2), 2, any)
      
      if(any(idx)) each.cl[idx, 'match.thresh'] <- 1
      return(each.cl)
    }
    
    if (nrow(each.cl) > 0) each.cl <- checkMatches(each.cl)

    
    
    name.match <- each.cl[,2] == 1
    if(any(name.match)){
      matched.cl <- matchCl(each.cl[name.match, ,drop=FALSE])
    } 
    
    # All cell lines that don't match GNE
    if(any(!name.match)) { 
      nonmatch.cl <- each.cl[!name.match, ,drop=FALSE]
      # Identify all cell line IDs that don't match
      rle.x <- rle(gsub("^.+?_", "", rownames(nonmatch.cl)))
      match.idx <- c()
      for(ref.idx in seq_along(rle.x$values)){
        # Individually prompt user to validate that the cell lines don't match
        ui.r <- NULL
        if(interactive){
          while(is.null(ui.r)){
            ui <- readline(paste0("GNE: ", colnames(nonmatch.cl)[1], 
                                  " / Ref: ", rle.x$values[ref.idx], 
                                  " - Does it match [Yes/No]: "))
            if(grepl("yes|y", ui, ignore.case = TRUE)) ui.r <- TRUE
            if(grepl("no|n", ui, ignore.case = TRUE)) ui.r <- FALSE 
          }
        } else {
          ui.r <- FALSE
        }
        
        
        # If the cell lines match, change the matching status to 1
        if(ui.r){
          cs2 <- cumsum(rle.x$length)
          cs1 <- c(1, cs2[-length(cs2)] + 1)
          mtmp <- matchCl(nonmatch.cl[cs1[ref.idx]:cs2[ref.idx],,drop=FALSE])
          
          if(exists("matched.cl")) matched.cl <- cbind(matched.cl, mtmp) else matched.cl <- mtmp
          match.idx <- c(match.idx, c(cs1[ref.idx]:cs2[ref.idx]))
        }
      }
      if(length(match.idx) > 0) nonmatch.cl <- nonmatch.cl[-match.idx, ,drop=FALSE]
    } else {
      nonmatch.cl <- matrix(nrow=1, ncol=1)
      rownames(nonmatch.cl) <- NULL
    }
    
    if(!exists("matched.cl")){
      matched.cl <- initializeConcordanceMatrix()
      matched.cl[1,] <- NA
    }
    
    list("match.anno"=rbind.fill(geno.mat, as.data.frame(matched.cl)),
         "nonmatch.anno"=rownames(nonmatch.cl))
  })
  
  
  match.anno.tmp <- do.call("rbind", 
                           lapply(all.geno.mats, function(x) x[['match.anno']][,1:4]))
  match.anno.df <- rbind(match.anno.df, match.anno.tmp)
  rownames(match.anno.df) <- names(top.conc)
  match.anno.df <- round(as.matrix(match.anno.df), 3)
  
  nonmatch.anno.tmp <- lapply(all.geno.mats, function(x) {
    boolean.nonmatch <- !is.null(x[['nonmatch.anno']])
    labels.nonmatch <- paste(x[['nonmatch.anno']], collapse="|")
    data.frame("stat"=as.logical(boolean.nonmatch), 
               "ids"=as.character(labels.nonmatch))
  })
  nonmatch.anno.tmp <- do.call("rbind", nonmatch.anno.tmp)
  nonmatch.anno.df <- rbind(nonmatch.anno.df, nonmatch.anno.tmp)
  rm(top.conc)
}

save(match.anno.df, nonmatch.anno.df, file="nAraw.GNE_matchdf.Rdata")
