#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 8                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o snpDens.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e snpDens.%N.%j.err.log      # File to which STDERR will be written

######################################
##
## This script takes the vcf containing a SUBSET of the 3039 samples
## and all positions with good sequencing quality (DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted)
## These files can be found across the multiple folders: 
## /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_*samples
## within pangloss. 
## input files: 
## 1. full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted.SNmono.BiSNPs.finalSamples.vcf.gz
## 2. ORF.Sace.coordinates.allChroms - file should be generated with preparingORFregions.sh
## 
## output files: 
## 1. allsites.25samples.noORF.noMiss.SNVs.demoHist.vcf.gz/vcfNoORFsSNVs.vcf
## 2. allsites.25samples.noORF.noMiss.500.SNVs.snpden.chr{1..16}
## the (1) is the clean vcf : no missing data and only SNVs
## the (2) is the list of 500bp (TO BE SPECIFIED) regions with the number of SNVs per window
## (2) is the input of XXXXXXX.sh which randomly selects 1SNP/window
## across windows 1000bp apart (TO BE SPECIFIED)
##
##
## by Isabel Alves - June 2024
##
######################################

start=`date +%s`

inputDir="/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples"
ORFDir="/shared/home/ialves/3039Sace/03-data"
# the vcf containing all the positions but only the 27 samples: 4 taiw, 3 sake, 7 mantou
# 7 wine and 4 dairy. 
inputVCF="full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.JubinExtracted.SNmono.BiSNPs.finalSamples.vcf.gz"

# checking if there is input dir and input vcf
if [ ! -f "$inputVCF" ] | [ ! -d "$inputDir" ]; 
then
echo "ERROR: no inputfile or no input dir. "
exit 0 
fi

cd $inputDir
outputVCForf="allsites.25samples.noORF.demoHist.vcf.gz"
# checking if there is ORF file 
if [ ! -f "$ORFDir/ORF.Sace.coordinates.allChroms" ] 
then
echo "ERROR: no ORF file present in $ORFDir/ORF.Sace.coordinates.allChroms."
exit 0 
fi

prefix=`echo $outputVCForf | sed 's/\(.*\).noORF.demoHist.vcf.gz/\1/'`
windowSize=500

module load vcftools
module load bcftools

echo "Computing total number of sites within: $inputVCF"
bcftools +counts $inputVCF

bcftools view $inputVCF -T ^$ORFDir/ORF.Sace.coordinates.allChroms -Oz -o $outputVCForf --threads 8

# checking if there is outputVCForf
if [ ! -f "$outputVCForf" ] 
then
echo "ERROR: VCF without coding sites (within ORFs) does not exist."
exit 0 
fi

echo "Computing total number of sites after removing the ORF: $outputVCForf"
bcftools +counts $outputVCForf
echo "Excluding missing data from: $outputVCForf"
bcftools view -i 'F_MISSING=0' $outputVCForf -Oz -o $prefix.noORF.noMiss.demoHist.vcf.gz --threads 8
noMissVCF="$prefix.noORF.noMiss.demoHist.vcf.gz"

# checking if there is noMissVCF
if [ ! -f "$noMissVCF" ] 
then
echo "ERROR: VCF without missing data does not exist."
exit 0 
fi

echo "Computing total number of sites after removing sites with missing data: $noMissVCF"
bcftools +counts $noMissVCF
vcftools --gzvcf $noMissVCF --SNPdensity $windowSize --out $prefix.noORF.noMiss.$windowSize

# checking if there is the file with genomic windows.
if [ ! -f "$prefix.noORF.noMiss.$windowSize.snpden" ] 
then
echo "ERROR: File with genomic windows does not exist."
exit 0 
fi

sed '1d' $prefix.noORF.noMiss.$windowSize.snpden > noheader.$prefix.noORF.noMiss.$windowSize.snpden

less $noMissVCF | grep -v "#" > vcfNoORFs.vcf 

# below we open the two files noheader.$prefix.$windowSize.snpden vcfNoORFs.vcf
# with the first on memory we go over all the positions within the vcf
# and checks in which interval in the file in memory each position falls and add +1
# this gives the counts of the number of positions across 500bp windows.
# window size is define above. These windows are obtained using vcftools
awk 'FNR == NR { chr[FNR] = $1; start[FNR] = $2; end[FNR] = $2+500; intCounts[FNR]=0; count++; next } 
FNR != NR {
    #print $0
    for (i=1; i<=count; i++)

        if($1==chr[i] && $2>=start[i] && $2<end[i]) {
            intCounts[i]++;
            #print chr[i],start[i],end[i],intCounts[i];
            break;
        };          
}
END {
    for (i=1; i<=count; i++)
        print chr[i],start[i],end[i],intCounts[i] 
}' noheader.$prefix.noORF.noMiss.$windowSize.snpden vcfNoORFs.vcf > $prefix.noORF.noMiss.$windowSize.siteden

# checking if there is the file with genomic windows.
if [ ! -f "$prefix.noORF.noMiss.$windowSize.siteden" ] 
then
echo "ERROR: File with site counts within genomic windows does not exist."
exit 0 
fi

echo "Removing monomorphic sites from: $noMissVCF"
outPrefix=`echo $noMissVCF | sed 's/\(.*\).demoHist.vcf.gz/\1/'`
bcftools view --min-ac 1 $noMissVCF -Oz -o $outPrefix.SNVs.demoHist.vcf.gz --threads 8
# checking if there is the vcf file only with SNVs
if [ ! -f "$outPrefix.SNVs.demoHist.vcf.gz" ] 
then
echo "ERROR: No SNV vcf."
exit 0 
fi

echo "Computing the number of SNVs within: $outPrefix.SNVs.demoHist.vcf.gz"
bcftools +counts $outPrefix.SNVs.demoHist.vcf.gz

less $outPrefix.SNVs.demoHist.vcf.gz | grep -v "#" > vcfNoORFsSNVs.vcf 
# below we open the two files noheader.$prefix.noORF.noMiss.$windowSize.snpden vcfNoORFs.vcf
# with the first on memory we go over all the positions within the vcf
# and checks in which interval in the file in memory each position falls and add +1
# this gives the counts of the number of positions across 500bp windows.
# window size is define above. These windows are obtained using vcftools
awk 'FNR == NR { chr[FNR] = $1; start[FNR] = $2; end[FNR] = $2+500; intCounts[FNR]=0; count++; next } 
FNR != NR {
    #print $0
    for (i=1; i<=count; i++)

        if($1==chr[i] && $2>=start[i] && $2<end[i]) {
            intCounts[i]++;
            #print chr[i],start[i],end[i],intCounts[i];
            break;
        };          
}
END {
    for (i=1; i<=count; i++)
        print chr[i],start[i],end[i],intCounts[i] 
}' noheader.$prefix.noORF.noMiss.$windowSize.snpden vcfNoORFsSNVs.vcf > $outPrefix.$windowSize.SNVs.snpden

# checking if there is the file containing the nb of SNVs counts per genomic window
if [ ! -f "$outPrefix.$windowSize.SNVs.snpden" ] 
then
echo "ERROR: File with SNV counts within genomic windows does not exist."
exit 0 
fi

for chr in {1..16}; 
do 
    grep chromosome${chr}$' ' $outPrefix.$windowSize.SNVs.snpden > $outPrefix.$windowSize.SNVs.snpden.chr$chr; 
done

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
