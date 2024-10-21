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
##  sbatch convertREFtoAA_5.sh -g /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples/allChrom.1SNP.25samples.noORF.noMiss.SNVs.demoHist \
##  -s samples25.txt -l list_of_strains_per_Clade.txt -t "onlyNonCoding" -d /shared/home/ialves/demoHist_yeast3039/03-data
##
##
############

start=`date +%s`
conda activate rForDemoInf
cd ~/demoHist_yeast3039/02-scripts
chmod +x convertREFtoAA_5.R 

while getopts "g:s:l:t:d:" option
do 
    case "${option}"
        in
        g)gtmatrix=${OPTARG};;
        s)samples=${OPTARG};;
        l)strlist=${OPTARG};;
        t)snptype=${OPTARG};;
        d)pathToData=${OPTARG};;
        
    esac
done

echo "Matrix of Positions: $gtmatrix";
echo "List Samples: $samples";
echo "List Strains: $strlist";
echo "snptype: $snptype";
echo "Path do data: $pathToData"

pathToPosM="$gtmatrix.chrpos"
pathToGTM="$gtmatrix.coded.GT"
pathToSNPdf="$pathToData/SNPs_$snptype.csv"
pathToCladeInfo="$pathToData/operationalTable_Full3039Sace_Clades_Ploidy_Aneuploidy.csv"
tmpDir=`dirname $pathToPosM`
pathToSamples="$tmpDir/$samples"
pathToStrains="$tmpDir/$strlist"

#----------------------
# Checking if files exist
#----------------------
FCOUNT=0
if [ ! -f $pathToPosM ];
then
    echo "ERROR: Position matrix doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $pathToGTM ];
then
    echo "ERROR: GT matrix doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $pathToSNPdf ]
then
    echo "ERROR: List of SNPs doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $pathToSamples ]
then
    echo "ERROR: List of samples doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $pathToCladeInfo ]
then
    echo "ERROR: Operational table doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
fi

if [ ! -f $pathToStrains ]
then
    echo "ERROR: Strains df doesn't exist."
else
    FCOUNT=$((FCOUNT+1))
fi


if [ $FCOUNT -eq 6 ]
then
    echo "Running convertREFtoAA.R on $pathToGTM and $pathToPosM matrices"
    echo "The list of SNPs is the followin: $pathToSNPdf"
    Rscript convertREFtoAA.R $gtmatrix $samples $strlist $snptype $pathToData
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
