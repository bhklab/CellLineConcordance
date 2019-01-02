library(scales)
PDIR <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/rdataFiles/cnv/output/GNE_corrPlots'
nA <- file.path(PDIR, "input", "nAraw.GNE_matchdf.Rdata")
nB <- file.path(PDIR, "input", "nBraw.GNE_matchdf.Rdata")

data.type <- 'cnv' # either cnv or snp
low.match.threshold <- c("cnv"=0.6,
                         "snp"=0.8)
max.num.low.match <- 20
match.col <- 'brown'
coi.col <- "red"  # Colour for points matched to cell lines of interest
coi <- NA

addAnnoLines <- function(anno.points, max.height=0.75, x.pos=1, 
                         dist.y=0.05, dir='right', segcol=alpha("black", 0.5)){
  anno.points <- anno.points[order(anno.points, decreasing = TRUE)]
  anno.names <- names(anno.points)
  if(dir == 'right'){
    x.pos.1 <- (x.pos + 1)
    x.pos.2 <- (x.pos + 2)
    x.text <- (x.pos + 2.1)
  } else if (dir == 'left'){
    x.pos.1 <- (x.pos - 1)
    x.pos.2 <- (x.pos - 2)
    x.text <- (x.pos - 2.1)
  } else {
    stop("Requires a 'right' or 'left' direction")
  }
  
  # Diagonal line from points
  segments(x0=rep(x.pos, length(anno.points)),
           y0=as.numeric(anno.points),
           x1=rep(x.pos.1, length(anno.points)),
           y1=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y),
           col=segcol)
  # Horizontal line from diagonal to name
  segments(x0=rep(x.pos.1, length(anno.points)),
           y0=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y),
           x1=rep(x.pos.2, length(anno.points)),
           y1=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y),
           col=segcol)
  # Name of Segments
  text(x=rep(x.text, length(anno.points)),
       y=seq(max.height, (max.height - dist.y*(length(anno.points)-1)), by=-dist.y),
       labels=anno.names, 
       adj=if(dir == 'right') 0 else 1)
}
plotDensity <- function(){
  if(mat.count == 1){ par(mar=c(0, 0, 4.1, 4.1)) } else { par(mar=c(0, 4.1, 4.1, 0)) }
  #par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  # Plot Shaded density plots for the match and mismatch distributions
  plot(if(mat.count == 1){ dens.match$y } else { ( dens.match$y * -1 ) }, 
       dens.match$x, 
       type="l", col=match.col,
       ylim=c(if(data.type=='cnv') -1 else 0,1),
       xlab="",
       ylab=if(data.type=='cnv' & mat.count == 2) 'Weighted Pearson Correlation' else if(data.type=='snp' & mat.count == 1) 'Concordance',
       axes=FALSE)
  
  # Add in coloured shades
  with(dens.match, polygon(x=if(mat.count == 1){ dens.match$y } else { ( dens.match$y * -1 ) }, 
                           y= dens.match$x, 
                           col=alpha(match.col, 0.5)))
  
  axis(side = 3, at=c(-100, 100), labels=FALSE) #Top  
  axis(side = 1, at=c(-100, 100), labels=FALSE) #Bottom
  if(mat.count == 1){
    axis(side = 4, at=c(-100, 100), labels=FALSE) #solid black bar to the right
  }
  if(mat.count == 2){
    axis(side = 2, at=c(-100, -1, -0.5, 0, 0.5, 1.0, 100), labels=TRUE) # solid black bar + labels to the left
  }
  
  abline(h = low.match.threshold[data.type], col="grey", lty=5)
  abline(h = quant.match.cutoff, col="dimgrey", lty=3)
  
  # Plot the individual points for match and mismatch points - highlighting the mismatch ones
  par(mar=c(0, 0, 4.1, 0))
  plot(rep(if(mat.count == 1) 1 else 9, length(unique(match.v))), 
       unique(match.v), 
       col=alpha(match.col, 0.2), 
       ylim=c(if(data.type=='cnv') -1 else 0,1), xlim=c(0,10), las=2,
       axes=FALSE,
       xlab='', ylab="",
       main=paste("hscr-a", if(mat.count == 1) "2" else "1", sep=""))
  #paste(gsub(".filename(.+)?", "", fn.headers, perl=TRUE), sep="", collapse=" - "))
  text(x=if(mat.count == 1) 3.1 else 6.9, 
       y=1.0, labels = paste("n = ", length(match.v), sep=""), 
       adj=if(mat.count == 1) 0 else 1)
  
  axis(side = 3, at=c(-100, 100), labels=FALSE) #Top  
  axis(side = 1, at=c(-100, 100), labels=FALSE) #Bottom
  if(mat.count == 1){
    axis(side = 2, at=c(-100, 100), labels=FALSE, col="grey") #gray bar to the left
  }
  if(mat.count == 2){
    axis(side = 4, at=c(-100, 100), labels=FALSE, col="grey") # gray bar to the right
  }
  
  
  # Add in low-match concordance "matching" points
  low.match <- match.v[which(match.v < cutoff.val)]
  if(length(low.match) < max.num.low.match) max.lowm <- length(low.match) else max.lowm <- max.num.low.match
  low.match <- sort(low.match)[c(0:max.lowm)]
  print(low.match)
  #if(!all(is.na(coi))) low.match <- match.v[which(names(match.v) %in% coi)]
  
  points(rep(if(mat.count == 1) 1 else 9, length(low.match)),  low.match, 
         pch=if(!all(is.na(coi))) 16 else 1,
         col=if(!all(is.na(coi))) coi.col  else match.col)
  if(length(low.match) > 0){
    addAnnoLines(anno.points=low.match, 
                 max.height=if(data.type=='cnv') (cutoff.val - 0.05) else 0.75,
                 x.pos=if(mat.count == 1) 1.15 else 8.95,
                 dir=if(mat.count == 1){  'right' } else { 'left' })
  }
  abline(h = low.match.threshold[data.type], col="grey", lty=5)
  abline(h = quant.match.cutoff, col="dimgrey", lty=3)
  
  if(mat.count == 2){
    par(mar=c(0,0,0,0))
    plot(0, type='n', xlim=c(0,10), ylim=c(0,10), axes=FALSE, ylab='', xlab='')
    legend(0.6,
           10, # places a legend at the appropriate place 
           c('matching cell lines'), # puts text in the legend
           lty=c(1), # gives the legend appropriate symbols (lines)
           lwd=c(2.5),col=c(match.col)) # gives the legend lines the correct color and width
  }
}

load(nA); nA <- match.anno.df
load(nB); nB <- match.anno.df
dens <- list("A"=nA, "B"=nB)


for(ds in colnames(nA)){
  print(paste0("Analyzing ", ds, "..."))
  pdf(file.path(PDIR, "output", paste0(ds, "_densityPlots.pdf")))
  layout(matrix(c(rep(c(0,3,4,4,2,2,1,0), 9),
                  c(0,0,0,0,5,5,0,0)),
                ncol=8, nrow=10, byrow=TRUE))
  for(nX in names(dens)){
    mat.count <- grep(nX, names(dens))
    match.v <- na.omit(dens[[nX]][,grep(ds, colnames(nA))])
    dens.match <- density(match.v, na.rm=TRUE) 
    quant.match.cutoff <- median(match.v) - (3*mad(match.v, center=median(match.v), constant=1.4826, na.rm=TRUE))
    if(quant.match.cutoff < low.match.threshold[data.type]) {cutoff.val <- low.match.threshold[data.type] } else { cutoff.val <- quant.match.cutoff}
    
    plotDensity()
  }
  dev.off()
}



