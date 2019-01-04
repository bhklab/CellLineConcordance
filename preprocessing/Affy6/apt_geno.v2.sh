#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#

/mnt/work1/users/bhklab/Tools/apt-1.16.1-x86_64-intel-linux/bin/apt-probeset-genotype \
-c /mnt/work1/users/bhklab/Tools/apt-1.16.1-x86_64-intel-linux/reference/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.cdf \
-a birdseed-v2 \
--read-models-birdseed /mnt/work1/users/bhklab/Tools/apt-1.16.1-x86_64-intel-linux/reference/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.birdseed.models \
--special-snps /mnt/work1/users/bhklab/Tools/apt-1.16.1-x86_64-intel-linux/reference/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.specialSNPs \
--out-dir /mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2015_09_GDSC/apt/apt_genotype_v2 \
--set-gender-method cn-probe-chrXY-ratio \
--chrX-probes /mnt/work1/users/pughlab/references/SNP_6.0/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.chrXprobes \
--chrY-probes /mnt/work1/users/pughlab/references/SNP_6.0/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.chrYprobes \
--write-models \
--cel-files /mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2015_09_GDSC/apt/filt.cel_list.txt
