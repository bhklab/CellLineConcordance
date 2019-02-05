#### Libraries ####
library(plyr)
library(scales)
library(data.table)

#### Functions ####
{
  getIdx <- function(idx.type, val, refdf, altdf){
    if(idx.type=='row'){
      grep(paste0("^", rownames(refdf)[val], "$"), rownames(altdf))
    } else if(idx.type == 'col'){
      grep(paste0("^", colnames(refdf)[val], "$"), colnames(altdf))
    }
  }

  # Functions
  cnvThresh <- function(x, cnv.threshold=0.8){
    x < cnv.threshold & x > -1
  }

  mergeCnv <- function(cnv.a, cnv.b){
    cnv.a <- unlist(strsplit(cnv.a, split = ","))
    cnv.b <- unlist(strsplit(cnv.b, split = ","))

    cnv.total <- unique(sort(c(cnv.a, cnv.b)))
    cnv.total <- paste(cnv.total, collapse=",")
    return(cnv.total)
  }

  matchCnvSnp <- function(cnv.a, snp.a){
    cnv.a <- unlist(strsplit(cnv.a, split = ","))
    snp.a <- unlist(strsplit(snp.a, split = ","))

    i <- intersect(cnv.a, snp.a)
    cnv.set <- setdiff(cnv.a, snp.a)
    snp.set <- setdiff(snp.a, cnv.a)

    return(list("cnv"=c("intersect"=length(i),
                        "cnv_only"=length(cnv.set)),
                "snp"=c("intersect"=length(i),
                        "snp_only"=length(snp.set)),
                "cnv_ids"=sort(cnv.set),
                "snp_ids"=sort(snp.set),
                "intersect"=sort(i)))
  }

  getVals <- function(cnv.snp.cnt, type='cnv'){
    sapply(cnv.snp.cnt, function(i){
      sapply(i, function(j){
        j[[type]]
      })
    })
  }

  mergeDsCnv <- function(data.str, cols, rows){
    cols <- unlist(cols)
    rows <- unlist(rows)

    data.str[data.str == ''] <- NA
    val <- na.omit(unique(sort(unlist(data.str[rows, cols]))))
    val <- paste(val, collapse=",")
    return(val)
  }
}

#### Main ####
options( stringsAsFactors = FALSE) ;
fn.headers <- c('cgp.filename', 'ccle.filename', "pfizer.filename.x", "pfizer2.filename.y", 'pfizer3.filename', 'gdsc.filename.x', 'gdsc2.filename.y')
names(fn.headers) <- c("cgp", "ccle", "pfizer", "pfizer2", "pfizer3", "gdsc", "gdsc2")

# Identify cell lines that are only found in one dataset
anno.dir <- "/Users/rquevedo/git/cnv_fingerprint/v2.0/data"
pdir <- '/Users/rquevedo/git/cnv_fingerprint/v2.0/summarize'
outdir <- file.path(pdir, "plots")
snp.dir <- "/Users/rquevedo/git/snp_fingerprint/environments"
cnv.dir <- '~/Onedrive/bhk_lab/results/phenotypes/data'

cnv.dir <- '/data/cnv'
outdir <- '/results/cnv'
snp.dir <- '/results/snp'  #codeocean
anno.dir <- "/data/ref"  #codeocean
setwd(outdir)
load(file.path(anno.dir,"merged_annotations.Rdata"))
multiple.cl.idx <- apply(cell.line.anno[,fn.headers], 1, function(x) length(which(is.na(x))) >= (length(fn.headers) - 1))

## Setting analysis flags
run.snp.analysis <- TRUE  # Generates a table of low-matching cell line identities between datasets
run.cnv.analysis <- TRUE # Generates a table of low-matching CNv concordances between datasets
visualize <- TRUE # Create stacked barplot visualizations

#### SNP Analysis ####
if(run.snp.analysis) {
  ### Load in match SNP data
  # Creates a dataframe containing all "Concordant GT and Discordant GT between similarly annotated cell lines"
  # [1] "y.df"              "low.match.list"    "cell.line.anno"    "matrix.perc.match"
  load(file.path(snp.dir, "plotPairwiseSnp.v2.Rdata"))
  load(file.path(anno.dir,"merged_annotations.Rdata"))

  # Tallies the number of low-matches for each dataset against every other dataset
  low.match.cnt <- lapply(low.match.list, function(x){
    low.cl <- apply(x, 1, function(each.row) which(each.row != 1 & !is.na(each.row)))
    low.cl.vals <- factor(unlist(lapply(low.cl, function(i) colnames(x)[i])))
    low.cl.cnt <- sapply(names(fn.headers), function(each.ds){
      ds.vals <- names(low.cl.vals)[as.character(low.cl.vals) == each.ds]
      tryCatch({paste(ds.vals, collapse=",")}, error=function(e){ NA })
    })

    #low.cl.cnt <- table(unlist(low.cl))
    #names(low.cl.cnt) <- colnames(x)[as.integer(names(low.cl.cnt))]
    return(data.frame(t(low.cl.cnt)))
  })

  low.match.cnt <- rbind.fill(low.match.cnt[names(low.match.cnt[[1]])])
  low.match.cnt[(low.match.cnt == 0)] <- NA
  rownames(low.match.cnt) <- colnames(low.match.cnt)

  low.match.cnt$length <- sapply(low.match.list[rownames(low.match.cnt)], nrow)
  low.match.cnt$names <- sapply(low.match.list[rownames(low.match.cnt)],
                                function(i) paste(rownames(i), collapse=", "))
}

#### CNV Analysis ####
if(run.cnv.analysis){
  ### Load in match CNV data
  # Create a list that contains only dataset-specific "Discordances" between matching annotations
  cnv.threshold <- 0.80
  datset.grep <-c("gdsc", "ccle", "cgp", "pfizer")
  datset.grep <- c("cgp", "ccle", "gdsc", "pfizer", "pfizer2", "pfizer3", "gdsc2")
  cnv.disc.match <- list()

  hscr.cnv.cnt <- lapply(c("hscra1", "hscra2"), function(each.hscr) {
    load(file.path(cnv.dir, paste0(each.hscr, ".match_filt.Rdata")))  #threshold.match.filt
    mismatch.cnv <- sapply(threshold.match.filt, function(x) any(x < cnv.threshold & x > -1))
    threshold.match.filt <- threshold.match.filt[mismatch.cnv]

    # For each dataset-allele, finds any cell lines that are discordant for matching annotations
    all.ds.cnv.list <- lapply(datset.grep, function(each.ds){
      ds.stat <- sapply(threshold.match.filt, function(x){
        match.idx <- grep(paste0(each.ds, "_"), rownames(x))
        cnv.ds.stat <- FALSE
        if(length(match.idx) > 0){
          cnv.ds.stat <- any(x[match.idx,] < cnv.threshold & x[match.idx,] > -1)
        }
        return(cnv.ds.stat)
      })
      return(threshold.match.filt[which(ds.stat)])
    })
    names(all.ds.cnv.list) <- datset.grep
    cnv.disc.match[[each.hscr]] <- all.ds.cnv.list

    all.ds.cnt <- lapply(names(all.ds.cnv.list), function(ref.ds.name) {
      each.ds.mm <- all.ds.cnv.list[[ref.ds.name]]

      # Get the count for each mismatched cell line for reference dataset 1 (ref.ds.name)
      ds.cl.cnt <- sapply(names(each.ds.mm), function(ref.cl.name){
        each.cl.mm <- each.ds.mm[[ref.cl.name]]
        row.idx <- grep(paste0("^", ref.ds.name, "_"), rownames(each.cl.mm))
        # for each comparison Dataset, determines if there's a mismatch value
        cl.cnt <- sapply(names(fn.headers), function(x, row.idx) {
          tryCatch({
            ds.idx <- grep(paste0("^", x, "_"), colnames(each.cl.mm))
            if(cnvThresh(each.cl.mm[row.idx, ds.idx])) ref.cl.name else 0  # Need to grab the cell line name here
            }, error=function(e){ 0 }) # If dataset does not have a value, return 0
          }, row.idx=row.idx)
        return(cl.cnt)
        })

      # Take a sum of all the mismatach to attribute mismatch values per dataset
      sum.ds.cl.cnt <- c("cgp" = NA)
      if(length(ds.cl.cnt) > 0) sum.ds.cl.cnt <- apply(ds.cl.cnt, 1, function(x){
        x <- x[(x != '0')]
        tryCatch({paste(x, collapse=",")}, error=function(e){ NA })
      })

      return(as.data.frame(t(sum.ds.cl.cnt)))
      })

    names(all.ds.cnt) <- names(all.ds.cnv.list)
    all.ds.cnt <- rbind.fill(all.ds.cnt[colnames(all.ds.cnt[[1]])])
    all.ds.cnt[(all.ds.cnt == 0)] <- NA
    rownames(all.ds.cnt) <- colnames(all.ds.cnt)

    all.ds.cnt$length <- sapply(all.ds.cnv.list[rownames(all.ds.cnt)], length)
    all.ds.cnt$names <- sapply(all.ds.cnv.list[rownames(all.ds.cnt)],
                               function(i) paste(names(i), collapse=", "))

    return(all.ds.cnt)
  })

  names(hscr.cnv.cnt) <- c("hscra1", "hscra2")
}

#### Intersect SNP and CNV analysis ####
if(run.snp.analysis & run.cnv.analysis){
  merged.cnv.cnt <- hscr.cnv.cnt[[1]][,c(1:length(fn.headers))]
  cnv.snp.cnt <- list()
  cnv.snp.ids <- list()

  for(each.row in c(1:nrow(merged.cnv.cnt))){
    for(each.col in c(1:ncol(merged.cnv.cnt))){
      mval <- mergeCnv(hscr.cnv.cnt[[1]][getIdx('row', each.row, merged.cnv.cnt, hscr.cnv.cnt[[1]]),
                                         getIdx('col', each.col, merged.cnv.cnt, hscr.cnv.cnt[[1]])],
                       hscr.cnv.cnt[[2]][getIdx('row', each.row, merged.cnv.cnt, hscr.cnv.cnt[[2]]),
                                         getIdx('col', each.col, merged.cnv.cnt, hscr.cnv.cnt[[2]])])
      merged.cnv.cnt[each.row, each.col] <- mval

      ids  <- c(rownames(merged.cnv.cnt)[each.row], colnames(merged.cnv.cnt)[each.col])
      ids <- sort(ids)
      cnv.snp.val <- matchCnvSnp(cnv.a=merged.cnv.cnt[each.row, each.col],
                                 snp.a=low.match.cnt[getIdx('row', each.row, merged.cnv.cnt, low.match.cnt),
                                                     getIdx('col', each.col, merged.cnv.cnt, low.match.cnt)])


      if(ids[1] != ids[2]) cnv.snp.cnt[[ids[1]]][[ids[2]]] <- cnv.snp.val[c('cnv', 'snp')]
      if(ids[1] != ids[2]) cnv.snp.ids[[ids[1]]][[ids[2]]] <- cnv.snp.val[c('cnv_ids', 'snp_ids', 'intersect')]
    }
  }
  cnv.snp.iddf <- lapply(cnv.snp.ids, function(x) do.call("rbind", x))

  fn.idx.headers <- c("gdsc", "cgp", "ccle", "pfizer")
  snp.idx.headers <- sapply(fn.idx.headers, function(i) grep(i, rownames(low.match.cnt)))
  fn.idx.headers <- sapply(fn.idx.headers, function(i) grep(i, colnames(merged.cnv.cnt)))
  summ.cnv.cnt <- matrix(ncol=length(fn.idx.headers), nrow=length(fn.idx.headers))
  colnames(summ.cnv.cnt) <- names(fn.idx.headers)
  rownames(summ.cnv.cnt) <- names(fn.idx.headers)
  summ.cnv.snp.val <- summ.cnv.cnt
  summ.cnv.snp.list <- list()

  for(each.col in c(1:ncol(summ.cnv.cnt))){
    for(each.row in c(1:nrow(summ.cnv.cnt))){
      colnm <- colnames(summ.cnv.cnt)[each.col]
      rownm <- rownames(summ.cnv.cnt)[each.row]
      print(paste(rownm, colnm, sep='-'))

      ds1.col <- grep(colnm, colnames(hscr.cnv.cnt[[1]]))
      ds1.row <- grep(rownm, rownames(hscr.cnv.cnt[[1]]))
      ds2.col <- grep(colnm, colnames(low.match.cnt))
      ds2.row <- grep(rownm, rownames(low.match.cnt))

      hscr.a1 <- mergeDsCnv(hscr.cnv.cnt[[1]], ds1.col, ds1.row)
      hscr.a2 <- mergeDsCnv(hscr.cnv.cnt[[2]], ds1.col, ds1.row)
      snp.data <- mergeDsCnv(low.match.cnt, ds2.col, ds2.row)

      summ.cnv.cnt[each.row, each.col] <- mergeCnv(hscr.a1, hscr.a2)

      ids  <- c(rownames(summ.cnv.cnt)[each.row], colnames(summ.cnv.cnt)[each.col])
      ids <- sort(ids)
      cnv.snp.val <- matchCnvSnp(cnv.a=summ.cnv.cnt[each.row, each.col],
                                 snp.a=snp.data)


      if(ids[1] != ids[2]) summ.cnv.snp.list[[ids[1]]][[ids[2]]] <- cnv.snp.val

    }
  }
}

#### Visualization ####
if(visualize){
  pdf(file.path(outdir, "karyo_geno_discPlot.pdf"), height = 5)
  cnv.snp.ds <- summ.cnv.snp.list
  #cnv.snp.ds <- cnv.snp.cnt

  cnv.x <- getVals(cnv.snp.ds, 'cnv')
  snp.x <- getVals(cnv.snp.ds, 'snp')
  box.cnv <- do.call("cbind", cnv.x)
  box.snp <- do.call("cbind", snp.x)

  id.x <- data.frame("ds1"=rep(names(cnv.snp.ds), sapply(cnv.snp.ds, length)),
                     "ds2"=colnames(box.snp))
  id.x <- id.x[c(nrow(id.x):1),]
  box.snp <- box.snp[,c(ncol(box.snp):1)]
  box.cnv <- box.cnv[,c(ncol(box.cnv):1)]

  dataset.pairs <- sum(sapply(cnv.snp.ds, length))

  split.screen(c(1,3))
  screen(1) #Identifier Screen
  par(mar=c(5.1, 4.1, 4.1, 0))
  cex.val <- 0.7
  plot(0, type='n', xlim=c(0,10), ylim=c(0, dataset.pairs), axes=FALSE, ylab='', xlab='')
  for(i in c(1:dataset.pairs)) text(x=1, y=i - 0.5, labels = id.x$ds1[i], cex=cex.val, adj = 0)
  for(i in c(1:dataset.pairs)) text(x=5, y=i - 0.5, labels = id.x$ds2[i], cex=cex.val, adj = 0)

  screen(2) #CNV Screen (Left, <0)
  par(mar=c(5.1, 0.5, 4.1, 0.25))
  bar_p <- barplot(box.cnv*-1, horiz=TRUE, xlim=c(-50, 0), ylab='', yaxt='n', xaxt='n')
  axis(side = 1, at=seq(0, -50, by=-10), labels=seq(0, 50, by=10), cex.axis=cex.val)
  axis(side = 3, tick=FALSE, at=-20, labels = "Karyotype", adj=0, cex.axis=cex.val)
  text(x = (-apply(box.cnv, 2, sum) - 1), y = bar_p + 0.25,
       labels=box.cnv[1,], adj=1, cex=cex.val)
  text(x = (-apply(box.cnv, 2, sum) - 1), y = bar_p - 0.25,
       labels=box.cnv[2,], adj=1, cex=cex.val)

  screen(3) #SNP Screen (Right, >0)
  par(mar=c(5.1, 0.25, 4.1, 0.5))
  bar_p <- barplot(box.snp, horiz=TRUE, xlim=c(0,50), ylab='', yaxt='n', cex.axis=cex.val)
  axis(side = 3, tick=FALSE, at=20, labels = "Genotype", adj=0, cex.axis=cex.val, cex.axis=cex.val)
  text(x = (apply(box.snp, 2, sum) + 1), y = bar_p + 0.25,
       labels=box.snp[1,], adj=0, cex=cex.val)
  text(x = (apply(box.snp, 2, sum) + 1), y = bar_p - 0.25,
       labels=box.snp[2,], adj=0, cex=cex.val)
  close.screen(all.screens=TRUE)
  dev.off()
}
# cnv.snp.iddf ## Details on intersect ids

# save(cnv.snp.iddf, file="intersectInfo/cnv_snp_idDF.Rdata")
cnp.snp.iddf <- lapply(names(cnv.snp.iddf), function(id_x){
  x <- cnv.snp.iddf[[id_x]]
  for(each_row in seq(nrow(x))){
    for(each_col in seq(ncol(x))){
      x[each_row, each_col] <- paste(unlist(x[each_row,each_col]), collapse=", ")

    }
  }
  x <- cbind(rep(id_x, nrow(x)), rownames(x), x)
  x
})
write.table(do.call("rbind", cnp.snp.iddf),
            # file="intersectInfo/cnv_snp_idDF.tsv",
            file=file.path(outdir, "cnv_snp_idDF.tsv"),
            sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)
