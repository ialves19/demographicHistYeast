#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o coRsfs.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e coRsfs.%N.%j.err.log      # File to which STDERR will be written

############
##
##
## this script should be run as:
## sbatch orrectingSFS_8.sh -i <inDir> -s <sitesFile> -m <modelFile>
##
## This script requires: 
##  1. an input Dir, 
##  2. model File,
##  3. sites file with the information about monomorphic, polymorphic and final dset. 
##
##  the input Dir: is where the sfs are. 
##  This script calls: correctingSFS_8.R
##  In the model file clade name should be in the right order of the model (.tpl)
##
##  Isabel Alves - March 2024
##
############


start=`date +%s`
. "/shared/home/ialves/anaconda3/etc/profile.d/conda.sh"
conda activate rForDemoInf
cd ~/demoHist_yeast3039/02-scripts
chmod +x correctingSFS_8.R

while getopts "i:s:m:" option
do 
    case "${option}"
        in
        
        i)indir=${OPTARG};;
        s)sites=${OPTARG};;
        m)model=${OPTARG};;
        
    esac
done

echo "Site information: $sites";
echo "Model file: $model";

#----------------------
# Checking if files exist
#----------------------
FCOUNT=0
if [ ! -f $indir/$model ]
then
    echo "ERROR: File with the model info doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
    modelTag=`head -1 $indir/$model`
fi

if [ $FCOUNT -eq 1 ] && [ -f $indir/${modelTag}_dSFS/$sites ]
then
    FCOUNT=$((FCOUNT+1))
else
    echo "ERROR: File with site info not found."
fi


if [ $FCOUNT -eq 2 ]
then
    echo "Running correctingSFS_8.R with input dir: $indir";
    echo "Correcting SFS for model: $model"
    echo "";
    Rscript correctingSFS_8.R $indir $sites $model
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