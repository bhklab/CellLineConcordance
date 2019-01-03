# CellLineConcordance

## cnv_analysis

### gCSI Analysis
  * **gCSI/createConcordanceMatrix.R:** Creates the base concordance matrices for genotypes between the Omni Array of gCSI with all the Affy6 array of CCLE/GDSC/Pfizer.  Saves the output as an Rdata structure for match.anno.df and nonmatch.anno.df
  * **gCSI/plotDrugLR.R:** Generates the heatmap of gCSI vs all other datasets: *genotype.conc.pdf*, as well as the copy-number concordance heatmap for nA and nB alleles: *[nBraw|nAraw].conc.pdf*

### CN Concordance Matrix
  * **cnvAvaFingerprinter.ignoretemp.R:** Used to generate the all-CL by all-CL CNV concordance matrix for A, B, and totaly copy ratio.
  * **cnvAvaFingerprinter.permutations.R:** Used to generate the all-CL by all-CL CNV concordance matrix for A, B, and totaly copy ratio. More focused on parallel processing and 

### Visualize CN of cell line pairs
  * **weightedNoisePlot.interactiveFunction:** Used to generate copy number figures for any two cell line pairs following quantile normalization.  Plots the total copy ratio track on the top panel, the euclidean distance between the two CN-profiles in the middle track, and the bottom track contains the variance for each segment.

### Visualize CN of cell line pairs
  * **plotPairwiseDatasetAnno.CNV.R:** Used to generate middle panel of Fig 4. Takes two datasets and generates a density distribution for matching cell line pairs (homonymous) and nonmatching cell lines (heteronomyous).  Annotates the cell lines that fall below a certain threshold.

### Summary metrics  
  * **summaryCclStats.discordant.R:** Used to generate Figure 2b.  Counts the number of homonymous cell line pairs that are discordant for their genotype, discordant for their karyotype, and discordant for both genotype and karyotype.
  * **summaryCclStats.concordant.R:** Used to generate Figure 2a.  Counts the number of heteronymous cell line pairs with concordant genotypes, then segregates them into proper match/mismtach category and intersects the list with ICLAC