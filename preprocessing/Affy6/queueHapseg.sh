#!/usr/bin/perl
use warnings;
use strict;
use Cwd;

my $dir = getcwd;
my $hapseg_out = '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2015_09_GDSC/hapseg/output';
my $quant_summ = '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2015_09_GDSC/apt/apt_summarize/quant-norm.pm-only.med-polish.expr.summary.txt';
my $bird_calls = '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2015_09_GDSC/apt/apt_genotype_v2/birdseed-v2.calls.txt';
my $bird_snp_models = '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/2015_09_GDSC/apt/apt_genotype_v2/birdseed-v2.snp-models.txt';


opendir(my $dir_fh, $dir) or die "Could not open directory, $dir, $!\n";

my $cnt = 1;
foreach my $file(sort(readdir($dir_fh))){
 
 if($file =~ /.*sub.*/){
   my $outfilename = $file;
   $outfilename =~ s/\_sub.+.txt/_$cnt\.sh/;
   $outfilename =~ s/^/hapseg\./;
   my $outfile = $dir . '/sh_scripts/' . $outfilename;
   open(my $out_fh, ">", $outfile) or die "Could not create file $outfilename, $!\n";
 
   my $filepath = $dir . '/' . $file;
   print $out_fh  <<BASH;
#!/bin/bash
#
#\$ -cwd
#\$ -S /bin/bash
#

module load R/3.0.2 

Rscript /mnt/work1/users/pughlab/src/HAPSEG/runHapSeg.hg19.R \\
$filepath \\
$hapseg_out \\
$quant_summ \\
$bird_calls \\
$bird_snp_models
BASH
   $cnt++;
 }
}

