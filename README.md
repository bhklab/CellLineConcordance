# CellLineConcordance

## phenotype_analysis
  * **gne_cnvDensityPlots.R:** Used to generate *Supplementary Figure 7: gCSI CNV-concordance density plots*.  Takes in the homonymous and heteronymous pairs of cell lines between gCSI and all other datasets and computes the matching and non-matching density plots for each pair of datasets.
  * **processPhenotype:** Used to generate a 3-panel phenotype display.  The top panel would indicate the Area Above the dose-response Curve (AAC) for all overlapping drugs between two cell lines. The lower panels would represent the comparison between expression values. 
  * **getAbcMatrices.R:** Computes the **A**rea **B**etween the dose-response **C**urves (ABC) for all overlaping drug concentrations between any two datasets.  Saves the resulting matrix as an RData file.
  * **gne_ccle.DrugPheno.R:** Used to generate *Figure 6: gCSI-CCLE CNV to drug phenotype comparison*. Takes the homonymous and non-homoynmous cell line pairs between gCSI and CCLE, separates into 3 exemplar cases per group, then plots the log-likelihood ratio for each cell line pair.

  * **cnvConc-drugAAC.R:** Computes the scatterplot relation betwee CNV concordance and delta-AAC values for each overlapping drug between GDSC and CCLE.