packages <- c("AnnotationDbi", "org.Hs.eg.db", "plyr",
              "vioplot", "PharmacoGx", "scales", "gtools",
              "intervals", "preprocessCore", "Hmisc",
              "ggplot2", "parallel")
tmp <- suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE))
args = commandArgs(trailingOnly=TRUE)

getMetaId <- function(pset, drug, cell){
  if(class(pset) == 'PharmacoSet'){
    ds.info <- pset@sensitivity$info
  } else if(class(pset) == 'data.frame'){
    ds.info <- pset
  } else {
    stop("Unknown data type: pset")
  }

  drug.idx <- grep(paste0("^", drug, "$"), ds.info$drugid, ignore.case = TRUE)
  cell.idx <- grep(paste0("^", cell, "$"), ds.info$cellid, ignore.case = TRUE)
  ds.id <- rownames(ds.info)[intersect(drug.idx, cell.idx)]
  ds.id
}

outdir="/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/abc/gcsi-ccle-abc/"
outdir='/data/pharmaco'

psets <- availablePSets()
ds.id1 <- psets$PSet.Name[grep(paste0("^", args[1], "$"), psets$PSet.Name, ignore.case = TRUE)]
ds.id2 <- psets$PSet.Name[grep(paste0("^", args[2], "$"), psets$PSet.Name, ignore.case = TRUE)]
#load("/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/abc/gcsi-ccle-abc/CCLE.Rdata")
gCSI <- downloadPSet(ds.id1)
#load("/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/abc/gcsi-ccle-abc/gCSI.Rdata")
CCLE <- downloadPSet(ds.id2)

dsA=gCSI
dsB=CCLE

dsA.info <- dsA@sensitivity$info
dsB.info <- dsB@sensitivity$info

cellids <- unique(sort(c(as.character(dsA.info$cellid),
                         as.character(dsB.info$cellid))))
drugids <- sort(c(unique(as.character(dsA.info$drugid)),
                  unique(as.character(dsB.info$drugid))))
drugids <- drugids[which(duplicated(drugids))]

dsA.drugs <- split(dsA.info, f=dsA.info$drugid)
dsB.drugs <- split(dsB.info, f=dsB.info$drugid)

abc.matrices <- lapply(drugids, function(drug.x){
  hills <-  list()
  s <- Sys.time()

  sapply(cellids, function(cell.i){
    # Gets the ID, conc, viability, and Hill Coeff for dataset A
    dsA.id <- getMetaId(pset=dsA.drugs[[drug.x]], drug=drug.x, cell=cell.i)[1]
    if(is.na(dsA.id)) dsA.id <- numeric()

    if(any(dsA.id == rownames(dsA@sensitivity$raw))){
      dsA.conc <- as.numeric(as.character(dsA@sensitivity$raw[dsA.id, , 'Dose']))
      dsA.viab <- as.numeric(as.character(dsA@sensitivity$raw[dsA.id, , 'Viability']))
      if(is.null(hills[[dsA.id]]) & !(all(is.na(dsA.conc)))){
        print(paste0("A: ", dsA.id, " - ", Sys.time(), " | ", round(Sys.time() - s, 3)))
        s <<- Sys.time()
        hills[[dsA.id]] <<- logLogisticRegression(conc=dsA.conc, viability=dsA.viab,
                                           conc_as_log = FALSE, viability_as_pct = TRUE,
                                           trunc = TRUE, verbose = TRUE)
      }
      dsA.hill <- hills[[dsA.id]]
    } else {
      dsA.id <- character()
      dsA.hill <- NULL
    }

    sapply(cellids, function(cell.j){
      # Gets the ID, conc, viability, and Hill Coeff for dataset B
      dsB.id <- getMetaId(pset=dsB.drugs[[drug.x]], drug=drug.x, cell=cell.j)[1]
      if(is.na(dsB.id)) dsB.id <- numeric()

      if(any(dsB.id == rownames(dsB@sensitivity$raw))){
        dsB.conc <- as.numeric(as.character(dsB@sensitivity$raw[dsB.id, , 'Dose']))
        dsB.viab <- as.numeric(as.character(dsB@sensitivity$raw[dsB.id, , 'Viability']))
        if(is.null(hills[[dsB.id]]) & !(all(is.na(dsB.conc)))){
          hills[[dsB.id]] <<- logLogisticRegression(conc=dsB.conc, viability=dsB.viab,
                                                   conc_as_log = FALSE, viability_as_pct = TRUE,
                                                   trunc = TRUE, verbose = TRUE)
        }
        dsB.hill <- hills[[dsB.id]]

      } else {
        dsB.id <- character()
        dsB.hill <- NULL
      }

      if((length(dsA.id) == 1) & (length(dsB.id) == 1) &
         (!is.null(dsA.hill)) & (!is.null(dsB.hill))){
        abc <- suppressWarnings(
          computeABC(conc1=dsA.conc, conc2 = dsB.conc,
                     Hill_fit1 = dsA.hill, Hill_fit2 = dsB.hill))
        round(abc/100, 3)
      } else {
        NA
      }
    })
  })
})

sapply(abc.matrices, class)

names(abc.matrices) <- drugids
attributes(abc.matrices)$colnames <- ds.id1
attributes(abc.matrices)$rownames <- ds.id2
outname <- paste0("abc.", tolower(ds.id1), "-",
                  tolower(ds.id2), ".Rdata")
save(abc.matrices, file=file.path(outdir, outname))


