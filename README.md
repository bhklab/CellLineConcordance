# CellLineConcordance

## cnv_analysis

  * **gCSI/createConcordanceMatrix.R:** Creates the base concordance matrices for genotypes between the Omni Array of gCSI with all the Affy6 array of CCLE/GDSC/Pfizer.  Saves the output as an Rdata structure for match.anno.df and nonmatch.anno.df
  * **gCSI/plotDrugLR.R:** Generates the heatmap of gCSI vs all other datasets: *genotype.conc.pdf*, as well as the copy-number concordance heatmap for nA and nB alleles: *[nBraw|nAraw].conc.pdf*