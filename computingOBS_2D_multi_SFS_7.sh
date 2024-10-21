#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o d2sfs.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e d2sfs.%N.%j.err.log      # File to which STDERR will be written

############
##
##
## this script should be run as:
## sbatch computingOBS_2D_multi_SFS_7.sh -i <inDir> -o <outDir> -s <samplesFile> -d <listofStrains> -m <modelFile> -g <GTm> 
## example: 
##          sbatch computingOBS_2D_multi_SFS_7.sh -i "/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples" \
##                  -o /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/fiveClades -s samples25.txt -d df_Clade_Strains.txt \
##                  -m model_TwoDomest_AsiOld_TaiwDip.txt -g <GTm>
##
## This script requires: 
##  1. an input Dir -i, 
##  2. an output Dir -o, 
##  3. a file with the samples harmonizedName -s,
##  4. a table with the samples/clade -d,
##  5. a genotype m already converted AA -g,
##  6. model File -m.
##
##  the input Dir: is where the output file of the convertREFtoAA.sh are. 
##  the output Dir: is where the SFS are gonna be generated. 
##  This script calls: computingOBS_2D_multi_SFS_7.R
##  the seed number is related to the randomly drawn allele at heterozygous
##  sites in the homozygous strains.
##  In case you want to treat the homozygous strains as heterozygous, 
##  change the MODEL FILE with the zygosity accordingly. 
##
##  Isabel Alves - March 2024
##
############


start=`date +%s`
. "/shared/home/ialves/anaconda3/etc/profile.d/conda.sh"
conda activate rForDemoInf
cd ~/demoHist_yeast3039/02-scripts
chmod +x computingOBS_2D_multi_SFS_7.R

while getopts "i:o:s:d:m:g:" option
do 
    case "${option}"
        in
        
        i)indir=${OPTARG};;
        o)outdir=${OPTARG};;
        s)samples=${OPTARG};;
        d)dfstrains=${OPTARG};;
        m)model=${OPTARG};;
        g)gtmatrix=${OPTARG};;
        
    esac
done

echo "Matrix of Positions: $gtmatrix";
echo "List Samples: $samples";
echo "Strain df: $dfstrains";
echo "Model tag: $model";

#----------------------
# Checking if files exist
#----------------------
FCOUNT=0
if [ ! -f $indir/$samples ];
then
    echo "ERROR: File with samples not found."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $indir/$dfstrains ];
then
    echo "ERROR: File with the strains/clade not found."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $indir/$gtmatrix ]
then
    echo "ERROR: genotype matrix doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $outdir/$model ]
then
    echo "ERROR: File with the model info doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ $FCOUNT -eq 4 ]
then
    echo "Running computingOBS_2D_multi_SFS_7.R with input dir: $indir and output dir: $outdir";
    echo "Computing SFS for model: $model"
    echo "";
    Rscript computingOBS_2D_multi_SFS_7.R $indir $outdir $dfstrains $gtmatrix $model $samples
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