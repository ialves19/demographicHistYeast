#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o demoII.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e demoII.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

wrkDir=~/demoHist_yeast3039/04-analysis
outputMatrices="demoInf_137Strains_SNVs"
cd $wrkDir

module load bcftools

# 1. getting the total number of monomorphic and variable sites after excluding missing data
echo "The total number of sites and the number of the variable sites is: "
bcftools view -i 'F_MISSING=0' full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted.SNmono.BiSNPs.vcf.gz | bcftools +counts

# 2. Excluding monomorphic sites 
echo "The total number of sites and the number of the variable sites is: "
bcftools view -i 'F_MISSING=0' full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted.SNmono.BiSNPs.vcf.gz | bcftools +fill-tags | bcftools view --min-ac 1 -Oz -o full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted.BiSNPs.vcf.gz --threads 8

# 3. Counting the total number of variable sites. 
echo "The number of variable sites is : "
bcftools +counts full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted.BiSNPs.vcf.gz

# 4. Exporting genotype and depth matrices 
bcftools query -f '[%GT\t]\n' full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted.BiSNPs.vcf.gz > ${outputMatrices}.GT
bcftools query -f '[%DP\t]\n' full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted.BiSNPs.vcf.gz > ${outputMatrices}.DP
bcftools query -f '%CHROM\t%POS\n' full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted.BiSNPs.vcf.gz > ${outputMatrices}.chrpos
#8.2 - Re-coding GT table 
cat ${outputMatrices}.GT | sed 's/|/\//g' | sed 's/\.\/\./-1/g' | sed 's/0\/0/0/g' | sed 's/0\/1/1/g' | sed 's/1\/0/1/g' | sed 's/1\/1/2/g' > ${outputMatrices}.coded.GT


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
