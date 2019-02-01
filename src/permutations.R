# Function: generateVarMatrix
# Purpose:  Gets all unique combinations for a given cell line boostrap and reports the 
#   squared difference for each bin
# Input:  A list element that contains all repeat permutations
# Returns:  Matrix: Contains all squared-difference for all bins (rows) and the corresponding std.dev weigh factor
generateVarMatrix <- function(uniq.cl.x=NULL, max.sd = 2, ...){
  require(utils)
  tcombn.mat <- NULL
  tcombn.adj.mat <- NULL
  row.sd <- NULL
  
  if(length(uniq.cl.x) > 1){
    uniq.combn <- combn(x = c(1:length(uniq.cl.x)), m = 2)
    
    #1)  Fills matrix M (m x n) with squared-diff values for :
    #   m: each genomic coordinate bin  (i.e. Chr1:5000-10000)
    #   n: matching cell lines from different subsets (i.e. sub10.sub3)
    tcombn.mat <- matrix(nrow=dim(total.rdata.list[[1]])[1], ncol=0)
    for(each.combn in c(1:dim(uniq.combn)[2])){
      id1 <- uniq.cl.x[uniq.combn[1,each.combn]]
      id2 <- uniq.cl.x[uniq.combn[2,each.combn]]
      combn.sd <- compareCn(total.rdata.list[[id1]][,seg.header], 
                            total.rdata.list[[id2]][,seg.header], 
                            median.norm=median.norm, quant.norm=quant.norm, ret='sd')
      
      # Append the squared-diff column to the matrix
      id1.n <- gsub("\\..+", "", id1)
      id2.n <- gsub("\\..+", "", id2)
      col.id <- paste(id1.n, id2.n, sep=".")
      tcombn.mat <- cbind(tcombn.mat, combn.sd)
      colnames(tcombn.mat)[dim(tcombn.mat)[2]] <- col.id
    }
    
    #2) Calculate the StdDev weighing factor:
    
    #     #a) 1/sd:
    #     row.sd <- apply(tcombn.mat, 1, function(x) 1/(sd(x) + 1))
    #     row.sd <- rescale(row.sd, to=c(0,1), 
    #                       from=c(1/(max.sd+1), 1) )  #(max.sd+1)
    
    #b) Squared 1/sd:
    row.sd <- apply(tcombn.mat, 1, function(x) 1/((sd(x) + 1)^2) )
    x <- rescale(row.sd, to=c(0,1), 
                 from=c(1/((max.sd +1)^2), 1 ))  #(max.sd+1)
    
    #     #c) Log Method:
    #     row.sd <- apply(tcombn.mat, 1, function(x) { 
    #       if(!is.na(sd(x))){ 
    #         log(1/(sd(x) + 1))
    #       } else { 
    #         log(1/(0 + 1))
    #       }})  #(sd(x) + 1)
    #     row.sd <- rescale(row.sd, to=c(0,1), 
    #                       from=c(log(1/(max.sd+1)), log(1/(0 + 1)) ))  #(max.sd+1)
    
    row.sd[which(row.sd < 0)] <- 0
    
    #3) Calculate weight-adjusted matrix:
    tcombn.adj.mat <- apply(tcombn.mat, 2, function(sdiff.col) return(sdiff.col * row.sd))
  }
  
  
  
  list("sdiff"=tcombn.mat, "adj.sdiff"=tcombn.adj.mat, "stdev"=row.sd)
}

