

removeDupMm <- function(x.df, col1, col2){
  # for each row in the Mismatch.df:
  dup.matches <- data.frame(matrix(ncol=2, nrow=0))
  for(each.row in 1:nrow(x.df)){
    match.vals <- as.character(x.df[each.row, c(col1, col2)])   #obtain the Mismatched cell lines
    match.rows <- which(apply(x.df, 1, function(x) length(which(match.vals %in% x)) > 1)) #Find the rows that contain that same mismatch
    if(!is.na(match.rows[2])){
      print(paste("Removing row: ", match.rows[c(2:length(match.rows))], sep=""))
      #rm.idx <- c(rm.idx, match.rows[c(2:length(match.rows))])
      dup.matches <- rbind(dup.matches, match.rows)
    }
  }
  
  # Order the row.indices dataframe so the lower row index is on the left
  dup.match.ord <- as.data.frame(t(apply(dup.matches, 1, function(x) if(x[1] > x[2]) c(x[2], x[1]) else x)))
  # Remove duplicates
  dup.match.ord <- dup.match.ord[order(dup.match.ord$X1),]
  dup.match.ord <- dup.match.ord[-which(duplicated(dup.match.ord)),]
  x.df <- x.df[-dup.match.ord[,2],]
  
  # Order the original mismatch.df so the "lowest" cell line is CL1
  x.df.ord <- as.data.frame(t(apply(x.df, 1, function(x) {
    if(x['cl1'] > x['cl2']) {
      c(x['X2'], x['X1'], x['value'], x['anno'], x['fingerprint'], x['cl2'], x['ds2'], x['cl1'], x['ds1'])
    } else x
  })))
  colnames(x.df.ord) <- colnames(x.df)
  x.df <- x.df.ord
  
  
  return(x.df)
}

getSelectedMismatchList <- function(input.list, all.ds.names, mismatch.m.df){
  x.list <- list()
  for(each.ds1 in all.ds.names){
    # Function to use a reference dataset and cell line (given by the excel supp.tables) and checks all other
    # datasets for the mismatched or unknown relation cell lines
    # Will return a vector of names with the Reference in Col1 and all misidentified in subsequent Cols
    x.list.temp <- apply(input.list, 1, function(x){
      v1.row.num <- grep(paste(each.ds1, x['V1'], sep="_"), mismatch.m.df$Var1)
      v2.row.num <- which(mismatch.m.df$Var2 %in% paste(all.ds.names, x['V2'], sep="_"))
      
      v1.v2.int <- intersect(v1.row.num, v2.row.num)
      
      return(c(mismatch.m.df[v1.v2.int[1], 'Var1'], mismatch.m.df[v1.v2.int, 'Var2']))
    })
    x.list[[each.ds1]] <- x.list.temp
  }
  return(x.list)
}

compressMismatchList <- function(input.list, all.ds.names, mismatch.m.df, mm.type, mismatch.dir){
  x.list <- list()
  for(each.ds in all.ds.names){
    # Compresses the list from the previous set into a Dataset-Reference dataframe and removes NA 
    # values where there is no cell line for the reference   
    x.list[[each.ds]] <- rbind.fill(lapply(input.list[[each.ds]], function(x){as.data.frame(t(x), stringsAsFactors=FALSE)}))
    if(!all(is.na(x.list[[each.ds]]))){
      colnames(x.list[[each.ds]]) <- c("Ref", 
                                       paste(rep("alt", (dim(x.list[[each.ds]])[2] - 1)),
                                             seq(1,(dim(x.list[[each.ds]])[2] - 1), by=1), sep=""))
    } else {
      colnames(x.list[[each.ds]]) <- c("Ref")
    }
    na.index <- which(is.na(x.list[[each.ds]]$Ref))
    x.list[[each.ds]] <- x.list[[each.ds]][-na.index,]
    
    #Print out a LaTeX file
    if(!is.null(dim(x.list[[each.ds]]))){
      bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
      print(xtable(x.list[[each.ds]]), include.rownames=FALSE,
            sanitize.colnames.function=bold,
            booktabs=T, file=file.path(mismatch.dir, mm.type, paste(each.ds, ".", mm.type, ".tex", sep="")))
    }
  }
  return(x.list)
}

