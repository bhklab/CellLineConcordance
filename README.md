# CellLineConcordance

## cnv_analysis

  * **gCSI/createConcordanceMatrix.R:** Creates the base concordance matrices for genotypes between the Omni Array of gCSI with all the Affy6 array of CCLE/GDSC/Pfizer.  Saves the output as an Rdata structure for match.anno.df and nonmatch.anno.df
  * **gCSI/plotDrugLR.R:** Generates the heatmap of gCSI vs all other datasets: *genotype.conc.pdf*, as well as the copy-number concordance heatmap for nA and nB alleles: *[nBraw|nAraw].conc.pdf*
  * **gCSI/plotDrugLR.R:** Generates the heatmap of gCSI vs all other datasets: *genotype.conc.pdf*, as well as the copy-number concordance heatmap for nA and nB alleles: *[nBraw|nAraw].conc.pdf*
  * **weightedNoisePlot.interactiveFunction:** Used to generate copy number figures for any two cell line pairs following quantile normalization.  Plots the total copy ratio track on the top panel, the euclidean distance between the two CN-profiles in the middle track, and the bottom track contains the variance for each segment.
  * **summaryCclStats.discordant.R:** Used to generate Figure 2b.  Counts the number of homonymous cell line pairs that are discordant for their genotype, discordant for their karyotype, and discordant for both genotype and karyotype.
  * **summaryCclStats.concordant.R:** Used to generate Figure 2a.  Counts the number of heteronymous cell line pairs with concordant genotypes, then segregates them into proper match/mismtach category and intersects the list with ICLAC