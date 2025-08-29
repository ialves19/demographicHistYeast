#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 8                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o bootnd.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e bootnd.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`
conda activate rForDemoInf
cd ~/demoHist_yeast3039/02-scripts
chmod +x convertREFtoAA_5.R 

module load bcftools

# 1. Defining file names 
wrkDir="/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples"
fName="allsites.25samples.noORF.noMiss.SNVs.demoHist.vcf.gz"
outName="allsites.25samples.noORF.noMiss.SNVs.demoHist"
# 2. Exporting genotype and depth matrices 
bcftools query -f '[%GT\t]\n' $wrkDir/$fName > $wrkDir/${outName}.GT

bcftools query -f '%CHROM\t%POS\n' $wrkDir/$fName > $wrkDir/${outName}.chrpos

cat $wrkDir/${outName}.GT | sed 's/|/\//g' | sed 's/\.\/\./-1/g' | sed 's/0\/0/0/g' | sed 's/0\/1/1/g' | sed 's/1\/0/1/g' | sed 's/1\/1/2/g' > $wrkDir/${outName}.coded.GT

if [ -f "$wrkDir/${outName}.coded.GT" ] && [ -f "$wrkDir/${outName}.chrpos" ];
then
    # 3. running script convertREFtoAA_5.sh
    sbatch convertREFtoAA_5.sh -g $wrkDir/${outName} -s samples25.txt -l list_of_strains_per_Clade.txt -t "all" -d /shared/home/ialves/demoHist_yeast3039/03-data
else
    echo "ERROR: no GT nor chrpos matrix found";
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