#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#

/mnt/work1/users/bhklab/Tools/apt-1.16.1-x86_64-intel-linux/bin/apt-probeset-summarize \
--cdf-file /mnt/work1/users/bhklab/Tools/apt-1.16.1-x86_64-intel-linux/reference/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.cdf \
--analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true \
--target-sketch /mnt/work1/users/bhklab/Tools/apt-1.16.1-x86_64-intel-linux/reference/gw6/lib/hapmap.quant-norm.normalization-target.txt \
--out-dir /mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2015_09_GDSC/apt/apt_summarize \
--cel-files /mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2015_09_GDSC/apt/filt.cel_list.txt
