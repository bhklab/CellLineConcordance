# CellLineConcordance

<<<<<<< HEAD
## snp_analysis
  * **gCSI/createMapping.R:** Creates a mapping between the gCSI provided metadata and the samplesheet that was used for preprocessing while all the ID's were still labelled as "Unk" for unknown.
  * **illuminaToAffyGT.R:** Converts the CRLMM genotype data (0/0, 0/1, 1/1) to the birdseed standard (0, 1, 2).  Overlaps the SNPs based on genomic location
  
  * **vcfToAffyGT.R:** Takes a standard VCF and converts it ot he Affy6 standard format for comparison of any sample to the Affymetrix 6 standard used for this analysis
  * **plotPairwiseDatasetAnno.R:** Creates a density plot of matching (homonymous) and non-matching (heteronymous) cell line pairs based on genotype concordance.
  * **subHeatmapPlotter.snp.R:** Plots individual sub-heatmaps of cell line pairs.  For each cell line, it'll create a heatmap of concordance to all isogenic lines and separate them into "Match" and "Mismatch" based on their annotations. Also genereates Rdata structure containing lists of the Match and Mismatch categories, as well as the complete matrix of concordance.
=======
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
>>>>>>> cnv_analysis
