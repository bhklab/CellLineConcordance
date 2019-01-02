annoGene <- function(genelist=NULL, mart){
  tryCatch({
    getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand'),
          filters=c('hgnc_symbol', 'ensembl_gene_id'),
          values=genelist,
          mart=mart)
  }, error=function(e){
    getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand'),
          filters=c('hgnc_symbol'),
          values=genelist,
          mart=mart)
  })
}

comparePhenotypes <- function(cl1, ds1, cl2, ds2, gen.plots=FALSE,
                              outdir=NA, loading.factors=20, expr.ids=NULL, input.type=NULL,
                              debug.mode=FALSE, auc.discordance=NA, mark.genes=NA, 
                              out.type='pdf', ...){
  if(debug.mode){
    cat(paste(paste('cl1 <- "', cl1, '"', sep=""),
              paste('ds1 <- "', ds1, '"', sep=""),
              paste('cl2 <- "', cl2, '"', sep=""),
              paste('ds2 <- "', ds2, '"\n', sep=""), sep="\n"))
  }
  
  uid <- toupper(paste(paste(ds1, cl1, sep="_"), paste(ds2, cl2, sep="_"), sep="-"))
  auc.mat <- getAucMat(cl1, ds1, cl2, ds2)
  expr.mat <- getExprMat(cl1, ds1, cl2, ds2)
  if(!is.null(expr.ids)) expr.mat <- subsetExprMat(expr.mat=expr.mat, 
                                                   expr.ids=expr.ids, 
                                                   input.type=input.type, ...) #expr.idx=expr.idx, rm.na=FALSE
  
  if(!all(is.na(auc.discordance)) & nrow(auc.mat) > 0){
    auc.disc <- as.matrix(sapply(auc.discordance, function(x){
      if(x %in% rownames(auc.mat)){
        disc.val <- abs(auc.mat[x,2] - auc.mat[x,1])  
      } else {
        disc.val <- NA
      }
      
      return(disc.val)
    }))
  } else {
    auc.disc <- NA
  }
  
  if((class(auc.mat) == 'matrix') | (class(expr.mat) == 'matrix')){
    if(gen.plots){
      if(out.type == 'pdf'){
        pdf(file.path(outdir, 
                      paste(paste(ds1, gsub("\\/", "-", cl1), sep="_"), 
                            paste(ds2, gsub("\\/", "-", cl2), sep="_"), 
                            "pdf", sep=".")),
            width=10)
      } else {
        png(file.path(outdir, 
                      paste(paste(ds1, gsub("\\/", "-", cl1), sep="_"), 
                            paste(ds2, gsub("\\/", "-", cl2), sep="_"), 
                            "png", sep=".")),
            height = 7, width=10, units = 'in', res = 300)
      }
      
    }
    split.screen(c(2,1))
    screen(1)
    if(dim(auc.mat)[1] > 1){
      barplot(t(auc.mat), beside=TRUE, las=2, ylim=c(0,1), 
              ylab="AUC", col=c(ds.col[ds1], ds.col[ds2]), cex.names=0.75)
      legend("topright", fill = c(ds.col[ds1], ds.col[ds2]), 
             border = FALSE, legend = c(ds1, ds2))
    } else {
      print(paste("Error: no overlapping drugs found between ", 
                  ds1, "-", cl1, " and ", 
                  ds2, "-", cl2, sep=""))
    }
    
    screen(2)
    lower.split <- split.screen(c(1,2))
    # Plot Expression PC2 Loading factors
    # screen(lower.split[1])
    # par(mar=c(5.1, 4.1, 0, 2.1))
    # if(dim(auc.mat)[1] > 1){
    #   plot(auc.mat, pch=19, col=alpha("grey", 0.9), xlim=c(0,1), ylim=c(0,1),
    #        xlab=paste(ds1, "AUC", sep=" "),
    #        ylab=paste(ds2, "AUC", sep=" "),
    #        sub=paste("Pearson Correlation: ", round(cor(auc.mat)[1,2],2), sep=""))
    # }
    
    
    
    
    if(dim(expr.mat)[1] > 1){
      if(!is.na(loading.factors)){
        expr.hugo <- getPc2LoadingFactors(expr.mat, top.x=loading.factors)
        row.idx <- sapply(names(expr.hugo), function(x) grep(x, rownames(expr.mat)))
      } else {
        row.idx <- match(paste(mark.genes, "at", sep="_"), rownames(expr.mat))
        row.idx <- row.idx[!is.na(row.idx)]
      }
      expr.cor <- round(cor(expr.mat)[1,2],2)
      
      
      if(ds1 == ds2) {
        ds1 <- 1
        ds2 <- 2
      }
      screen(lower.split[2])
      par(mar=c(5.1, 4.1, 0, 1))
      plot(expr.mat, pch=19, col=alpha("grey", 0.3),
           xlab=paste(ds1, "expression", sep=" "),
           ylab=paste(ds2, "expression", sep=" "),
           xlim=c(0,15), ylim=c(0,15), las=1,
           sub=paste("Pearson Correlation: ", expr.cor, sep=""))
      points(x=expr.mat[row.idx, ds1],
             y=expr.mat[row.idx, ds2],
             col="red", pch=16)
      lines(x = c(0,100), y = c(0,100), lty=2)
      
      
      # Code to create KS-Plots in base plots obtained from:
      # https://rpubs.com/mharris/KSplot
      sq.residuals <- getSqResiduals(expr.mat, analysis='diagonal')
      
      rnd.smooth.spline <- smooth.spline(x=round(sq.residuals$xval,1),
                                         y=round(sq.residuals$residual,1))
      rnd.smooth.spline$y[which(rnd.smooth.spline$y < 0)] <- 1
      sample1 <- rep(rnd.smooth.spline$x, ceiling(rnd.smooth.spline$y))
      sample2 <- rnd.smooth.spline$x
      cdf1 <- ecdf(sample1)
      cdf2 <- ecdf(sample2)
      minMax <- seq(min(sample1, sample2), max(sample1, sample2), length.out=length(sample1)) 
      x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
      y0 <- cdf1(x0) 
      y1 <- cdf2(x0) 
      
      ks.pval <- ks.test(x=sample1, 
                         y=sample2, alternative="two.sided")$p.val
      screen(lower.split[1])
      par(mar=c(5.1, 4.1, 0, 2.1))
      plot(0, type="n", xlim=c(0,20), ylim=c(0,(7^2)),
           xlab="Expression", ylab="Squared residuals")
      text(x = 16, y = 48, labels = paste("p.val:", round(ks.pval,3)), cex=0.7, col="black")
      # segments(x0 = sq.residuals$xval, y0 = rep(0, length(sq.residuals$residual)), 
      #          x1 = sq.residuals$xval, y1 = sq.residuals$residual, col=alpha("black", 0.25))
      points(x = sq.residuals$xval, y=sq.residuals$residual, col=alpha("black", 0.25))
      lines(rnd.smooth.spline, col="red")
      subplot(plotCdf(cdf1, cdf2, x0, y0, y1), 
              x = 16, y=38, size = c(0.75, 0.75))
      
    } else {
      expr.cor <- NA
      print(paste("Error: no overlapping expression set found between ", 
                  ds1, "-", cl1, " and ", 
                  ds2, "-", cl2, sep=""))
    }
    close.screen(all.screens=TRUE)
    if(gen.plots) dev.off()
  }
  
  return(list("auc"=auc.disc,
              "expr"=expr.cor,
              "deltaAuc"=sort(abs(apply(auc.mat, 1, diff))),
              "loadingFactors"=expr.hugo))
}

getSqResiduals <- function(expr.mat, analysis='x'){
  if(analysis == 'x'){
    r.val <- (expr.mat[,2] - expr.mat[,1])^2
    sq.residuals <- data.frame("raw.x" = expr.mat[,1],
                               "raw.y" = expr.mat[,2],
                               "xval" = NA,
                               "yval" = NA,
                               "residual" = r.val)
  } else if (analysis == 'y'){
    r.val <- (expr.mat[,1] - expr.mat[,2])^2
    sq.residuals <- data.frame("raw.x" = expr.mat[,1],
                               "raw.y" = expr.mat[,2],
                               "xval" = NA,
                               "yval" = NA,
                               "residual" = r.val)
  } else if(analysis == 'diagonal'){
    b.val <- expr.mat[,2] - ((-1)*expr.mat[,1])
    x.intersect <- b.val/2
    y.intersect <- x.intersect
    
    r.val <- (expr.mat[,2] - y.intersect)^2 + (expr.mat[,1] - x.intersect)^2
    direction.val <- rep(1, length(r.val))
    direction.val[which(expr.mat[,2] < y.intersect)] <- -1
    r.val <- r.val * direction.val
    
    sq.residuals <- data.frame("raw.x" = expr.mat[,1],
                               "raw.y" = expr.mat[,2],
                               "xval" = x.intersect,
                               "yval" = y.intersect,
                               "residual" = r.val)
  } else {
    cat(paste("analysis type specified:", analysis, "\n"))
    cat("Argument parameters for analysis:
        \tanalysis = 'x'   # Squared difference in reference to X
        \tanalysis = 'y'   # Squared difference in reference to Y
        \tanalysis = 'diagonal'   # Squared difference in reference to diagonal of f(x)=x")
  }
  return(sq.residuals)
}

getPc2LoadingFactors <- function(expr.mat, top.x=40){
  print("PC")
  print(paste0("Dim: ", dim(expr.mat)))
  print(expr.mat[1:5,1:5])
  pca.expr <- prcomp(expr.mat)
  pca.expr$x <- pca.expr$x[order(pca.expr$x[,"PC2"]),]
  expr.len <- nrow(pca.expr$x)
  pca.expr$x <- pca.expr$x[c(1:top.x,(expr.len - top.x):expr.len),]
  expr.hugo <-  mapIds(org.Hs.eg.db,
                       keys=gsub("_at", "", rownames(pca.expr$x)),
                       keytype="ENSEMBL",
                       column="SYMBOL",
                       multiVals="first")
  return(expr.hugo)
}

getIdx <- function(cl, ds){
  cl.idx <- grep(paste("^", cl, "$", sep=""), 
                 colnames(ds), ignore.case=TRUE)
  if(length(cl.idx) == 0){
    print(paste("Error: Could not find ", cl, " in dataset", sep=""))
  }
  
  return(cl.idx)
}

plotCdf <- function(cdf1, cdf2, x0, y0, y1){
  plot(cdf1, verticals=TRUE, do.points=FALSE, col="red",
       xlab='', ylab='', main='',  yaxt='n', axes=FALSE)
  plot(cdf2, verticals=TRUE, do.points=FALSE, col="grey", add=TRUE)
  segments(x0, y0, x0, y1, col="red", lty="dotted") 
  axis(side = 2, at = c(0, 0.5, 1), labels=c(0, 0.5, 1), las=2, cex.axis=0.75)
  axis(side = 1, at = seq(0,16, by=4), labels=seq(0,16, by=4), las=1, cex.axis=0.75)
}

getAucMat <- function(cl1, ds1, cl2, ds2){
  cl1.idx <- getIdx(cl1, drug.auc.list[[ds1]])
  cl2.idx <- getIdx(cl2, drug.auc.list[[ds2]])
  
  merged.cl.mat <- matrix(nrow=0, ncol=1)
  if(all(sapply(list(cl1.idx, cl2.idx), function(x) length(x) > 0))){
    merged.cl.mat <- t(smartbind(drug.auc.list[[ds1]][,cl1.idx],
                                 drug.auc.list[[ds2]][,cl2.idx]))
    colnames(merged.cl.mat) <- c(ds1, ds2)
    na.rows <- apply(merged.cl.mat, 1, function(x) !any(is.na(x)))
    merged.cl.mat <- merged.cl.mat[which(na.rows),]
  }
  
  return(merged.cl.mat)
}

getExprMat <- function(cl1, ds1, cl2, ds2){
  cl1.idx <- getIdx(cl1, expr.list[[ds1]])
  cl2.idx <- getIdx(cl2, expr.list[[ds2]])
  
  merged.expr.mat <- matrix(nrow=0, ncol=1)
  if(all(sapply(list(cl1.idx, cl2.idx), function(x) length(x) > 0))){
    expr1.mat <- exprs(expr.list[[ds1]])[,cl1.idx, drop=FALSE]
    expr2.mat <- exprs(expr.list[[ds2]])[,cl2.idx, drop=FALSE]
    merged.expr.mat <- cbind(expr1.mat, expr2.mat[,1][match(rownames(expr1.mat), 
                                                            rownames(expr2.mat))])
    colnames(merged.expr.mat) <- c(ds1, ds2)
    na.rows <- apply(merged.expr.mat, 1, function(x) !any(is.na(x)))
    merged.expr.mat <- merged.expr.mat[which(na.rows),]
  }
  
  return(merged.expr.mat)
}

subsetExprMat <- function(expr.mat, expr.ids, input.type='ENTREZID'){
  conv.ens.id <- mapIds(org.Hs.eg.db,
                        keys=as.character(expr.ids),
                        column="ENSEMBL",
                        keytype=input.type,
                        multiVals="first")
  expr.idx <- match(paste(conv.ens.id, "at", sep="_"), rownames(expr.mat))
  expr.mat <- expr.mat[expr.idx,]
  
  na.idx <- which(apply(expr.mat, 1, function(x) all(is.na(x))))
  if(length(na.idx) > 0) expr.mat <- expr.mat[-na.idx,]
  return(expr.mat)
}
