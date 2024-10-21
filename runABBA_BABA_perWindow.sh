#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o ABBApw.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e ABBApw.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

# loading packages
module load bcftools
conda activate r
# Defining working dir and input files
wrkDir=/shared/home/ialves/demoHist_yeast3039/04-analysis/ABBA-BABA
cd $wrkDir
mkdir ABBA_BABA_perWindow
cd ABBA_BABA_perWindow
inDir=/shared/home/ialves/3039Sace/03-data
inFile="full1279_diploidNoAneu_biSNPs_MQ30_MQRankSum-12.5_excHet1_updated.vcf.gz"
outputPrefix="AW_AD_ED.mac1.Dsuite"

# file with 1279 strain names and respective clades
inSamples="1279strainClades.txt"


cp ../cladeList.txt .
#cat $inSamples | sort -k 2 -n > 1279strainClades.ordered.txt
cat cladeList.txt | while read line; do grep "$line" 1279strainClades.txt; done > strains_of_cladeList.txt
cat strains_of_cladeList.txt | awk '{print $1}' > samplesToKeep.tmp
bcftools view -S samplesToKeep.tmp ${inDir}/${inFile} | bcftools +fill-tags | bcftools view --min-ac 1 -Oz -o $outputPrefix.vcf.gz --threads 8

bcftools +counts $outputPrefix.vcf.gz

#bcftools query -f '[%GT\t]\n' AW_AD_ED.mac1.Dsuite.vcf.gz > AW_AD_ED.mac1.Dsuite.GT
bcftools query -f '%CHROM\t%POS\n' $outputPrefix.vcf.gz > $outputPrefix.chrpos


Rscript 

outputRalt="$outputPrefix.allSNVs.alt"
outputRref="$outputPrefix.allSNVs.ref"

if [ ! -f "$outputRalt" ] || [ ! -f "$outputRref" ];
then

echo "ERROR: alt or ref files missing."
exit 0;
fi

bcftools index -t $outputPrefix.vcf.gz

bcftools view -R $outputRalt $outputPrefix.vcf.gz -Oz -o $outputPrefix.ALT.vcf.gz --threads 8 
bcftools view -R $outputRref $outputPrefix.vcf.gz -Oz -o $outputPrefix.REF.vcf.gz --threads 8 

gunzip -ck $outputPrefix.REF.vcf.gz | \
awk 'BEGIN {OFS="\t"} {if ( $1 ~ /^chromosome/ ) print $0,"0/0"; else if ( $1 ~ /^#CHROM/ ) print $0,"SPAR"; else print $0;}' > $outputPrefix.REF_SPAR.vcf

gunzip -ck $outputPrefix.ALT.vcf.gz | \
awk 'BEGIN {OFS="\t"} {if ( $1 ~ /^chromosome/ ) { t = $4; $4 = $5; $5 = t; print;} else { print $0; }}' > $outputPrefix.ALT_TMP.vcf

cat $outputPrefix.ALT_TMP.vcf | \
awk 'BEGIN {OFS="\t"} {if ( $1 ~ /^chromosome/ ) print $0,"0/0"; else if ( $1 ~ /^#CHROM/ ) print $0,"SPAR"; else print $0;}' > $outputPrefix.ALT_SPAR_TMP.vcf

cat $outputPrefix.ALT_SPAR_TMP.vcf | \
sed 's/|/\//g' | sed 's/1\/1/2\/2/g' | sed 's/0\/0/1\/1/g' | sed 's/2\/2/0\/0/g' > $outputPrefix.ALT_SPAR_AAconverted.vcf

bcftools view $outputPrefix.REF_SPAR.vcf -Oz -o $outputPrefix.REF_SPAR.vcf.gz
bcftools index -c $outputPrefix.REF_SPAR.vcf.gz
bcftools view $outputPrefix.ALT_SPAR_AAconverted.vcf -Oz -o $outputPrefix.ALT_SPAR_AAconverted.vcf.gz
bcftools index -c $outputPrefix.ALT_SPAR_AAconverted.vcf.gz
bcftools concat -a $outputPrefix.REF_SPAR.vcf.gz $outputPrefix.ALT_SPAR_AAconverted.vcf.gz -Oz -o $outputPrefix.SPAR_AAconverted.vcf.gz --threads 8 

$HOME/Dsuite/Build/Dsuite Dinvestigate -w 10,5 AW_AD_ED.mac1.Dsuite.SPAR_AAconverted.vcf.gz strains_of_cladeList.txt test_FD_WW.txt
$HOME/Dsuite/Build/Dsuite Dtrios AW_AD_ED.mac1.Dsuite.SPAR_AAconverted.vcf.gz strains_of_cladeList.txt

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