#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o d1sfs.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e d1sfs.%N.%j.err.log      # File to which STDERR will be written

############
##
##
## this script should be run as:
## sbatch computingOBS_1D_SFS_6.sh -d <wrkDir> -s <samplesFile> -f <listofStrains> -z <zygosityFile> -g <GTm> 
##
## This script requires: 
##  1. the genotype matrix, 
##  2. the sample file, 
##  3. the table with strains/clade, and
##  4. the table with the zygosity 
##
##  to be in the same dir: wrkDir
##  This script calls: computingOBS_1D_SFS_6.R
##  the seed number is related to the randomly drawn allele at heterozygous
##  sites in the homozygous strains.
##  In case you want to treat the homozygous strains as heterozygous, 
##  change the the table with the zygosity accordingly. 
##
##  Isabel Alves - March 2024
##
############


start=`date +%s`
. "/shared/home/ialves/anaconda3/etc/profile.d/conda.sh"
conda activate rForDemoInf
cd ~/demoHist_yeast3039/02-scripts
chmod +x computingOBS_1D_SFS_6.R

while getopts "d:s:f:z:g:" option
do 
    case "${option}"
        in
        
        d)wrkdir=${OPTARG};;
        s)samples=${OPTARG};;
        f)dfstrains=${OPTARG};;
        z)zig=${OPTARG};;
        g)gtmatrix=${OPTARG};;
        
    esac
done

echo "Matrix of Positions: $gtmatrix";
echo "List Samples: $samples";
echo "Strain df: $dfstrains";
echo "Zygosity df: $zig";

#----------------------
# Checking if files exist
#----------------------
FCOUNT=0
if [ ! -f $wrkdir/$samples ];
then
    echo "ERROR: File with samples not found."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $wrkdir/$dfstrains ];
then
    echo "ERROR: File with the strains/clade not found."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $wrkdir/$gtmatrix ]
then
    echo "ERROR: genotype matrix doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $wrkdir/$zig ]
then
    echo "ERROR: File with zygosity info doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ $FCOUNT -eq 4 ]
then
    echo "Running computingOBS_1D_SFS_6.R on $wrkdir and $gtmatrix matrices";
    echo "";
    Rscript computingOBS_1D_SFS_6.R $wrkdir $samples $dfstrains $zig $gtmatrix
else 
    echo "ERROR: no input files found." 
fi

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