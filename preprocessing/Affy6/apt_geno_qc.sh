#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#

/mnt/work1/users/bhklab/Tools/apt-1.16.1-x86_64-intel-linux/bin/apt-geno-qc \
--cdf-file /mnt/work1/users/bhklab/Tools/apt-1.16.1-x86_64-intel-linux/reference/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.cdf \
--qca-file /mnt/work1/users/pughlab/references/SNP_6.0/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.r2.qca \
--qcc-file /mnt/work1/users/pughlab/references/SNP_6.0/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.r2.qcc \
--chrX-probes /mnt/work1/users/pughlab/references/SNP_6.0/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.chrXprobes \
--chrY-probes /mnt/work1/users/pughlab/references/SNP_6.0/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.chrYprobes \
--cel-files /mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2015_09_GDSC/apt/cel_list.txt \
--out-file /mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2015_09_GDSC/apt/apt_geno_qc/results.txt

