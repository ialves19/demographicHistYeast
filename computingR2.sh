#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 8                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o rSquared.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e rSquared.%N.%j.err.log      # File to which STDERR will be written


######################################
##
## This script takes the vcf containing only biallelic SNPs outside ORFs:
## allsites.*samples.noORF.noMiss.SNVs.demoHist.vcf.gz
## within the following directory:
## /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_*samples
## in pangloss. 
## input files: 
## 1. allsites.*samples.noORF.noMiss.SNVs.demoHist.vcf.gz
## 
## output files: 
## 1. rsquared.allSNVs.geno.ld
## 
## the (1) is the output of vcftools with R2 values in windows of XX bp. no missing data and only SNVs
##
## Window size and output name needs to be SPECIFIED. 
## The scrpit should be run like: 
## sbatch computingR2.sh <windowsize> <input Name> <output Name>
##
## by Isabel Alves - Oct 2024
##
######################################


start=`date +%s`

# loading packages
module load vcftools

cd /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_137samples
# arguments to pass to VCFTOOLS
windSize=$1
inFile=$2
outFile=$3

# run vcftools
vcftools --gzvcf $inFile --geno-r2 -ld-window-bp $windSize --out $outFile


end=`date +%s`
runtime=$((end-start))
#days
D=$((runtime / 60 / 60 / 24))
# all hours
H=$((runtime / 60 / 60))
H=$((H % 24))
# all minutes
M=$((runtime / 60))
M=$((M % 60))
# all seconds
S=$((runtime % 60))
echo "This task took: $D days, $H hours, $M minutes and $S seconds."