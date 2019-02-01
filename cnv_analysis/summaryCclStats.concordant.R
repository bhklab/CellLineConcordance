library(plyr)
library(scales)

### -------------------- Functions
processTallies <- function(tally.df, col.list, split.ds){
  if(class(tally.df) == 'list') tally.df <- rbind.fill(tally.df) else tally.df <- t(tally.df)
  relevant.col.ids <- which(names(col.list) %in% colnames(tally.df))
  tally.df <- tally.df[,names(col.list)[relevant.col.ids], drop=FALSE]
  rownames(tally.df) <- names(split.ds)
  # tbl.cnts <- tbl.cnts[rev(rownames(tbl.cnts)),,drop=FALSE]
  tally.df[is.na(tally.df)] <- 0
  return(as.data.frame(tally.df))
}

#Assemble the NA/unformated table data into a structured matrix for stacked barplots
postProcessCnts <- function(all.tbl.cnts, all.ds.list){
  na.idx <- which(sapply(all.tbl.cnts, function(x) all(is.na(x))))
  names(all.tbl.cnts) <- names(all.ds.list)
  all.tbl.cnts <- all.tbl.cnts[-na.idx]
  rename.row.ids <- lapply(names(all.tbl.cnts), function(x) data.frame("a"=rep(x, nrow(all.tbl.cnts[[x]])),
                                                                       "b"=rownames(all.tbl.cnts[[x]])))
  all.tbl.cnts <- rbind.fill(all.tbl.cnts)
  all.tbl.cnts <- cbind(do.call("rbind", rename.row.ids), all.tbl.cnts)
  all.tble.cnts <- t(as.matrix(all.tbl.cnts))
  return(all.tble.cnts)
}



### -------------------- Pre-processing
# Load the recent output and data-structure from subHeatmapPlotter.snp.R that was reformatted on May 24-2017
#load("/Users/rquevedo/git/snp_fingerprint/environments/plotPairwiseSnp.v2.Rdata")    # same SNP matrix as subHeatmap.May24-2017.Rdata

# "matrix.perc.match" "mm.ds.list" "unknown.cl.list" "conc.list" "anno.mismatch.m.filt.df"
load("/Users/rquevedo/Desktop/bhk_lab/results/snp_fingerprinting/allSnps/subHeatmap.May24-2017.Rdata")  

options( stringsAsFactors = FALSE) ;

col.id <- list('match'="#fee5d9",
               'unknown'="#fcae91",
               'Possible match'="#fb6a4a",
               'mismatch'="#cb181d")
iclac.col <- list("0"="#91003f",
                  "1"="#e7298a",
                  "2"="#c994c7")


x123 <- anno.mismatch.m.filt.df[,c("ds1", "cl1", "ds2", "cl2", "iclac1", "iclac2", "concId")]
x123 <- as.data.frame(apply(x123, 2, as.character))
all.ds.names <- unique(sort(unlist(x123[,c("ds1", "ds2")])))
all.ds.list <- lapply(all.ds.names, function(each.ds){
  ds.idx <- apply(x123, 1, function(each.row) any(each.ds %in% each.row) )
  temp.x <- x123[which(ds.idx),]
  
  # Move the DS of interest to the left column (ds1, cl1) and the comparator to the right column (ds2, cl2)
  switch.idx <- which(!temp.x$ds1 %in% each.ds)
  temp.x.ds2df <- temp.x[switch.idx, c("ds2", "cl2", 'iclac2')]
  temp.x.ds1df <- temp.x[switch.idx, c("ds1", "cl1", 'iclac1')]
  
  temp.x[switch.idx, ]$cl1 <- as.character(temp.x.ds2df$cl2)
  temp.x[switch.idx, ]$ds1 <- as.character(temp.x.ds2df$ds2)
  temp.x[switch.idx, ]$iclac1 <- as.integer(temp.x.ds2df$iclac2)
  temp.x[switch.idx, ]$cl2 <- as.character(temp.x.ds1df$cl1)
  temp.x[switch.idx, ]$ds2 <- as.character(temp.x.ds1df$ds1)
  temp.x[switch.idx, ]$iclac2 <- as.integer(temp.x.ds1df$iclac1)
  
  return(temp.x)
})
names(all.ds.list) <- all.ds.names

### -------------------- Assemble stacked-barplot matrices of counts
iclac.anno.cnts <- lapply(names(all.ds.list), function(each.ds.df.id){
  each.ds.df <- all.ds.list[[each.ds.df.id]]
  split.ds <- split(each.ds.df, each.ds.df$ds2)
  
  # Remove repeating segments that were processed in previous iterations
  rm.ids <- names(all.ds.list)[c(0:(match(each.ds.df.id, names(all.ds.list)) - 1))]
  if(length(rm.ids) > 0) split.ds <- split.ds[-which(names(split.ds) %in% rm.ids)]

  iclac.tallies <- tryCatch({
    #Tally Iclac between each dataset
    iclac.tallies <- lapply(split.ds, function(ds.df) {
      ds.df <- ds.df[which(ds.df$concId == 'mismatch'),]  # Subset for only mismatches
      as.data.frame(t(as.matrix(table( apply(ds.df, 1, function(x) {
        sum(as.integer(x[c('iclac1', 'iclac2')]), na.rm=TRUE)
      })))))
    })
    iclac.tallies <- processTallies(iclac.tallies, iclac.col, split.ds)
    
    iclac.tallies
  }, error=function(e){NA})
  
  tbl.cnts <- tryCatch({
    #Tally the match/mismatch/possible matches between each dataset
    tbl.cnts <- sapply(split.ds, function(x) as.data.frame(t(as.matrix(table(x$concId))), drop=FALSE))
    tbl.cnts <- processTallies(tbl.cnts, col.id, split.ds)
    
    tbl.cnts
  }, error=function(e){NA})
  
  return(list("iclac"=iclac.tallies,
              "anno.tbl"=tbl.cnts))
})



all.anno.cnts <- postProcessCnts(lapply(iclac.anno.cnts, function(x) x[['anno.tbl']]), all.ds.list)
all.anno.cnts[is.na(all.anno.cnts)] <- 0
all.iclac.cnts <- postProcessCnts(lapply(iclac.anno.cnts, function(x) x[['iclac']]), all.ds.list)
all.iclac.cnts[is.na(all.iclac.cnts)] <- 0



### -------------------- Visualization
pdf("~/Desktop/annoCclConc.pdf")
split.screen(matrix(c(0, 0.7, 0, 1,
                      0.7, 1, 0, 1), byrow=TRUE, ncol=4))
screen(1)
par(mar=c(5.1, 8, 4.1, 0.5))
bar.vals <- barplot(all.anno.cnts[c(3:6),], horiz = TRUE, yaxt='n', 
                    col=as.character(col.id[rownames(all.anno.cnts)[3:6]]), las=1,
                    xlim=c(0, 80), cex.axis=0.7)
legend(x = 50, y = max(bar.vals) - 1, fill = unlist(col.id), 
       legend = names(col.id), cex=0.7)
axis(side = 2, at = bar.vals, labels = all.anno.cnts['b',], tick=FALSE, las=1, cex.axis=0.7, line=0)
axis(side = 2, at = bar.vals, labels = all.anno.cnts['a',], tick=FALSE, las=1, cex.axis=0.7, line=2.5)

screen(2)
par(mar=c(5.1, 0.5, 4.1, 2.1))
bar.vals <- barplot(all.iclac.cnts[c(3:5),], horiz = TRUE, yaxt='n', 
                    col=as.character(iclac.col[as.character(rownames(all.iclac.cnts)[3:5])]), 
                    las=1, xlim=c(0, 40), cex.axis=0.7)
legend(x = 15, y = max(bar.vals) - 1, fill = unlist(iclac.col), 
       legend = paste(names(iclac.col), "in ICLAC", sep=" "), cex=0.7)

close.screen(all.screens=TRUE)
dev.off()



### -------------------- Enrichment analysis
getContigencyTable <- function(x, col.id="all"){
  #x <- t(all.anno.cnts)[11,]
  cl.a <- x['a']
  cl.b <- x['b']
  
  if(col.id == 'all'){
    t.mismatch.cnt <- sum(as.integer(x[3:length(x)]))
  } else if(any(names(x) == col.id)){
    t.mismatch.cnt <- as.integer(x[match(col.id, names(x))])
  } else {
    stop(paste0("No matching column id: ", paste(names(x), collapse=", ")))
  }
  
  
  cl.a.idx <- grep(paste0("^", cl.a, ".filename"), colnames(cell.line.anno))
  cl.b.idx <- grep(paste0("^", cl.b, ".filename"), colnames(cell.line.anno))
  intersect.tbl <- cbind(cell.line.anno[,cl.a.idx],
                         cell.line.anno[,cl.b.idx])
  nonNa.cnt <- apply(intersect.tbl, 2, function(y) length(which(!is.na(y))))
  nonNa.cnt <- prod(nonNa.cnt)
  
  na.idx <- apply(intersect.tbl, 1, function(y) all(!is.na(y)))
  match.rm <- length(which(na.idx))
  nonNa.cnt <- nonNa.cnt - match.rm
  return(data.frame("match"=t.mismatch.cnt,
                    "mismatch"=(nonNa.cnt - t.mismatch.cnt)))
}


### Phased out because it would get us crucified
#########
# contigency.tbl <- do.call("rbind", apply(t(all.anno.cnts), 1, getContigencyTable, col.id="match"))
# rownames(contigency.tbl) <- paste(t(all.anno.cnts)[,1], t(all.anno.cnts)[,2], sep="-")
# contigency.tbl <- contigency.tbl[c(1:3, 7:8, 12),]
# Xsq.match <- chisq.test(contigency.tbl)
# 
# contigency.tbl <- do.call("rbind", apply(t(all.anno.cnts), 1, getContigencyTable, col.id="mismatch"))
# rownames(contigency.tbl) <- paste(t(all.anno.cnts)[,1], t(all.anno.cnts)[,2], sep="-")
# contigency.tbl <- contigency.tbl[c(1:3, 7:8, 12),]
# Xsq.mismatch <- chisq.test(contigency.tbl)
# 
# dim(Xsq.match$stdres)
# qnorm(1-((0.05/12)/2))
# 
# 
# Xsq.match$stdres
# Xsq.mismatch$stdres
# Xsq.mismatch$stdres / Xsq.match$stdres
