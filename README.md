# CellLineConcordance

## snp_analysis
  * **gCSI/createMapping.R:** Creates a mapping between the gCSI provided metadata and the samplesheet that was used for preprocessing while all the ID's were still labelled as "Unk" for unknown.
  * **illuminaToAffyGT.R:** Converts the CRLMM genotype data (0/0, 0/1, 1/1) to the birdseed standard (0, 1, 2).  Overlaps the SNPs based on genomic location
  
  * **vcfToAffyGT.R:** Takes a standard VCF and converts it ot he Affy6 standard format for comparison of any sample to the Affymetrix 6 standard used for this analysis
  * **plotPairwiseDatasetAnno.R:** Creates a density plot of matching (homonymous) and non-matching (heteronymous) cell line pairs based on genotype concordance.
  * **subHeatmapPlotter.snp.R:** Plots individual sub-heatmaps of cell line pairs.  For each cell line, it'll create a heatmap of concordance to all isogenic lines and separate them into "Match" and "Mismatch" based on their annotations. Also genereates Rdata structure containing lists of the Match and Mismatch categories, as well as the complete matrix of concordance.
  * **snpOneVsAllFingerprinter.R:** Takes a single sample and compares it against a reference panel (Affy6 Datasets).