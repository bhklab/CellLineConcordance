###################
#### Functions ####
plotTransition <- function(cols, ccle.spl, allele){
  plot(0, type='n', xlim=c(1,3), ylim=c(-1,1), xaxt='n', ylab='correlation (r)', las=1)
  summ.val <- lapply(names(cols), function(lvl){
    ccle.tmp <- ccle.spl[[lvl]]
    with(ccle.tmp, points(x=rep(1, nrow(ccle.tmp)), y=r_germline_CCLE_HC_vs_GDSC_WES, pch=16, col=cols[lvl]))
    with(ccle.tmp, points(x=rep(3, nrow(ccle.tmp)), y=r_somatic_CCLE_HC_vs_GDSC_WES, pch=16, col=cols[lvl]))
    
    if(allele == 'A1'){
      with(ccle.tmp, points(x=rep(2, nrow(ccle.tmp)), y=A1, pch=16, col=cols[lvl]))
      with(ccle.tmp, segments(x0 = rep(1, nrow(ccle.tmp)), y0 = r_germline_CCLE_HC_vs_GDSC_WES, 
                              x1 = rep(2, nrow(ccle.tmp)), y1 = A1, col=cols[lvl]))
      with(ccle.tmp, segments(x0 = rep(2, nrow(ccle.tmp)), y0 = A1, 
                              x1 = rep(3, nrow(ccle.tmp)), y1 = r_somatic_CCLE_HC_vs_GDSC_WES, col=cols[lvl]))
    } else if(allele=='A2') {
      with(ccle.tmp, points(x=rep(2, nrow(ccle.tmp)), y=A2, pch=16, col=cols[lvl]))
      with(ccle.tmp, segments(x0 = rep(1, nrow(ccle.tmp)), y0 = r_germline_CCLE_HC_vs_GDSC_WES, 
                              x1 = rep(2, nrow(ccle.tmp)), y1 = A2, col=cols[lvl]))
      with(ccle.tmp, segments(x0 = rep(2, nrow(ccle.tmp)), y0 = A2, 
                              x1 = rep(3, nrow(ccle.tmp)), y1 = r_somatic_CCLE_HC_vs_GDSC_WES, col=cols[lvl]))
    } else {
      stop("Not A1 or A2")
    }
    
    
    
    do.call(cbind, lapply(ccle.tmp[,c('r_germline_CCLE_HC_vs_GDSC_WES', 
                                      allele, 
                                      'r_somatic_CCLE_HC_vs_GDSC_WES')], summary))
  })
  legend("bottom", legend = c("Germline match", 
                              "Considerable drift (r.somatic < 0.75)",
                              "Substantial drift (r.somatic < 0.6)",
                              "Germline mismatch"), fill = cols)
  names(summ.val) <- names(cols)
  summ.val
}
plotVals <- function(vals, karyo="Karyotype (A2)"){
  plot(0, axes=FALSE, type='n', xlim=c(1,3), ylim=c(0.5,4.5), ylab='', xlab='')
  sapply(1:ncol(vals), function(x){
    sapply(1:nrow(vals), function(y){
      text(x, nrow(vals) - y + 1, labels=vals[y,x], cex=0.8)
    })
  })
  abline(v=c(1.5, 2.5))
  axis(side = 2, at=c(1:4), labels=rev(rownames(vals)), las=2, cex.axis=0.8, tick = FALSE)
  axis(side = 1, at = c(1,2,3), labels=c("Germline", karyo, "Somatic"), las=2)
}

###################
#### Load Data ####
ccldir <- '~/git/CellLineConcordance/cnv_analysis/data'
refdir <- '~/git/CellLineConcordance/ref'
outdir <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/cnv_analysis/ccle2_comparison'
setwd(outdir)

load(file.path(ccldir, 'hscra1.weightedCorrDf.Rdata')); hscra1 <- weighted.corr.df
load(file.path(ccldir, 'hscra2.weightedCorrDf.Rdata')); hscra2 <- weighted.corr.df
load(file.path(refdir, 'merged_annotations.Rdata'));
ccle2 <- read.table(file.path(refdir, 'ccle2_stable3.csv'), sep=",",
                    header=TRUE, stringsAsFactors = FALSE, check.names = FALSE);

#################
#### Map IDs ####
ccle2$ids <- sapply(ccle2$CCLE_ID, function(x) strsplit(x, split="_")[[1]][[1]])
mapped.ids <- t(sapply(ccle2$ids, function(x){
  idx <- grep(paste0("^", x, "$"), cell.line.anno$unique.cellid.atoz, ignore.case=TRUE)
  id <- cell.line.anno[idx, 'unique.cellid']
  if(x=='TT') {
    tt.idx <- which(id == 'TT')
    c(id[tt.idx], idx[tt.idx],
      cell.line.anno[idx[tt.idx], 'ccle.filename'],
      cell.line.anno[idx[tt.idx], 'gdsc.filename.x']) 
  } else if(length(id) == 0) {
     c(NA, NA, NA, NA)
  } else {
    c(id, idx, 
      cell.line.anno[idx, 'ccle.filename'],
      cell.line.anno[idx, 'gdsc.filename.x'])
  }
}))
mapped.ids <- as.data.frame(mapped.ids)
colnames(mapped.ids) <- c("mapped.ids", "mapped.idx",
                          "mapped.ccle", "mapped.gdsc")
mapped.ids$mapped.idx <- as.integer(as.character(mapped.ids$mapped.idx))
ccle2 <- cbind(ccle2, mapped.ids)

##############################
#### Map CNV Discordances ####
ccle2$mapped.gdsc <- gsub(".cel", "", as.character(ccle2$mapped.gdsc), ignore.case = TRUE)
ccle2$mapped.ccle <- gsub(".cel", "", as.character(ccle2$mapped.ccle), ignore.case = TRUE)
cnv.corr <- apply(ccle2, 1, function(cl.row){
  a1 <- hscra1[cl.row['mapped.ccle'], cl.row['mapped.gdsc']]
  if(is.null(a1)) a1 <- NA
  a2 <- hscra2[cl.row['mapped.ccle'], cl.row['mapped.gdsc']]
  if(is.null(a2)) a2 <- NA
  c(a1, a2)
})
cnv.corr <- as.data.frame(t(cnv.corr))
colnames(cnv.corr) <- c("A1", "A2")
ccle2 <- cbind(ccle2, cnv.corr)


#########################
#### Estimate counts ####
# Number of cell lines not assayed using SNP6
table(apply(ccle2, 1, function(x) all(is.na(c(x['mapped.ccle'], x['mapped.gdsc']))))) # 11 TRUE

# Number of cell lines not assayed in both datasets
table(apply(ccle2, 1, function(x) any(is.na(c(x['mapped.ccle'], x['mapped.gdsc']))))) # 24 TRUE

# Number of cell lines with viable karyotypes in both datasets
table(apply(ccle2, 1, function(x) any(is.na(c(x['A1'], x['A2']))))) # 129 TRUE; 524 FALSE

# Number of cell lines with high somatic (>0.75) and low karyotype (<0.8)
table(apply(ccle2, 1, function(x){
  any(c(x['A1'] < 0.8 & x['r_somatic_CCLE_HC_vs_GDSC_WES'] > 0.75,
        x['A2'] < 0.8 & x['r_somatic_CCLE_HC_vs_GDSC_WES'] > 0.75))
}))



snu <- c('SNU1', 'SNU16', 'SNU5')
idx <- unlist(sapply(snu, grep, cell.line.anno$unique.cellid))
cell.line.anno[which(cell.line.anno$unique.cellid == c('SNU1', 'SNU16', 'SNU5')),]

# Number of cell lines not found in both datasets
table(apply(ccle2, 1, function(x) any(is.na(c(x['mapped.ccle'], x['mapped.gdsc']))))) # 24 TRUE
table(apply(ccle2, 1, function(x) all(is.na(c(x['mapped.ccle'], x['mapped.gdsc']))))) # 24 TRUE
which(apply(ccle2, 1, function(x) all(is.na(c(x['mapped.ccle'], x['mapped.gdsc'])))))
# Number of cell lines missing A1 or A2
table(apply(ccle2, 1, function(x) any(is.na(c(x['A1'], x['A2']))))) # 129 TRUE

##################################
#### Plot based on categories ####
ccle.spl <- split(ccle2, f=ccle2$comments)
names(ccle.spl) <- c("match", "considerable", "mismatch", "substantial")

cols <- c("grey", "#fdcc8a", "#fc8d59", "#d7301f")
names(cols) <- c("match", "considerable", "substantial", "mismatch")

## Allele A
pdf("a1_ccle2-comp.pdf", height = 10, width = 12)
split.screen(c(1,2))
screen(1)
a1.scrn <- split.screen(matrix(c(0, 1, 0.3, 1.0,
                                 0, 1, 0, 0.3), byrow=TRUE, ncol=4))
screen(a1.scrn[1]); par(mar=c(0, 6.1, 4.1, 2.1))
summ.val <- plotTransition(cols, ccle.spl, 'A1')

screen(a1.scrn[2]); par(mar=c(10, 6.1, 0.5, 2.1))
vals <- t(sapply(summ.val, function(x) x['Median',]))
plotVals(vals, "Karyotype (A1)")

#close.screen(all.screens=TRUE)
#dev.off()

## Allele B
#pdf("a2_ccle2-comp.pdf", height = 10)
screen(2)
a2.scrn <- split.screen(matrix(c(0, 1, 0.3, 1.0,
                                 0, 1, 0, 0.3), byrow=TRUE, ncol=4))
screen(a2.scrn[1]); par(mar=c(0, 6.1, 4.1, 2.1))
summ.val <- plotTransition(cols, ccle.spl, 'A2')

screen(a2.scrn[2]); par(mar=c(10, 6.1, 0.5, 2.1))
vals <- t(sapply(summ.val, function(x) x['Median',]))
plotVals(vals, "Karyotype (A2)")

close.screen(all.screens=TRUE)
dev.off()
