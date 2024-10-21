#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o demoI.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e demoI.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

#inputDir="/mnt/viseg/Projects/FullSaceSNPmatrix"
wrkDir=~/demoHist_yeast3039/03-data
cd $wrkDir

module load bcftools

bcftools view -S samples31.txt full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.SNPs.vcf.gz | \
bcftools view -R SNPs_Dsuite.txt | bcftools +fill-tags | bcftools --min-ac 1 -Oz -o demographicInf_31samples_ANCeqREF_Dsuite.vcf.gz

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

