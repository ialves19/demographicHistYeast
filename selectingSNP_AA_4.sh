#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o rTest.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e rTest.%N.%j.err.log      # File to which STDERR will be written

############
##
##
## this script should be run as:
## sbatch selectingSNP_AA_4.sh <all or onlyNonCoding>
##
############

conda activate rForDemoInf
cd ~/demoHist_yeast3039/02-scripts
chmod +x selectingSNP_AA_4.R
# checking files 

aaFile=~/ancestralState/03-Data/full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.chromosome1.JubinExtracted.SNmono.BiSNPs.AA
annoFile=~/demoHist_yeast3039/03-data/annotation/full3039Matrix.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.ann_4SIFT_SIFTannotations.xls

FCOUNT=0
if [ ! -f $aaFile ];
then
echo "ERROR: Ancestral file doesn't exist."
else
FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $annoFile ];
then
echo "ERROR: Annotation file doesn't exist."
else
FCOUNT=$((FCOUNT+1))
fi

if [ $FCOUNT -eq 2 ]
then
	#export R_LIBS="/shared/home/ialves/R/x86_64-redhat-linux-gnu-library/4.3"
	Rscript selectingSNP_AA_4.R $1
else 
	echo "ERROR: no input files found." 
fi

