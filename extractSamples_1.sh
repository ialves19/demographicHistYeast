#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o demoI.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e demoI.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

inputDir="/shared/home/ialves/ancestralState/03-Data/allSites"
wrkDir="/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_137samples"
outDir="/shared/home/ialves/demoHist_yeast3039/04-analysis"
cd $wrkDir

module load bcftools

sed 1d list_of_strains_per_Clade.txt | cut -d$'\t' -f2 > strainsToKeep.tmp

bcftools concat ${inputDir}/full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.chromosome{1..16}.JubinExtracted.SNmono.BiSNPs.vcf.gz | \
bcftools view -S strainsToKeep.tmp | bcftools +fill-tags | bcftools view -Oz -o ${outDir}/full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted.SNmono.BiSNPs.vcf.gz --threads 8



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
