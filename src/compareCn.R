compareCn <- function(cn1, cn2, median.norm=0, quant.norm=0, ret='sd'){
  if(median.norm==1){
    print("Median normalizing the segments")
    avg.median <- mean(c(median(cn1), median(cn2)))
    
    cn1 <- ((cn1 - median(cn1)) + avg.median)
    cn2 <- ((cn2 - median(cn2)) + avg.median)
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