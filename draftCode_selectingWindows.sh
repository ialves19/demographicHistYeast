#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 8                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o 1SNVWind.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e 1SNVWind.%N.%j.err.log      # File to which STDERR will be written

######################################
##
## This script takes the vcf containing a SUBSET of the 3039 samples
## and ONLY SNVs without missing data and outside ORFs. 
##
## These files can be found across the multiple folders: 
## /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_*samples
## within pangloss. 
## input files: 
## 1. allsites.25samples.noORF.noMiss.SNVs.demoHist.vcf.gz
## 2. allsites.25samples.noORF.noMiss.500.SNVs.snpden.chr{1..16} - this is created with gettingSNPs_density.sh
## 
## output files: 
## 1. allChrom.1SNP.25samples.noORF.noMiss.SNVs.demoHist.vcf.gz
## 2. allChromosomes.1SNPperWindow.txt
## the (1) is the clean vcf : no missing data, only SNVs and ONLY ONE SNV per window
## the (2) is the list of genomic windows (500bp <- TO BE SPECIFIED) regions with the number of SNVs per window
## and the number of the SNVs randomly picked up. 
##
## by Isabel Alves - June 2024
##
######################################

start=`date +%s`

module load vcftools
module load bcftools

wrkDir=/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal
subDir="SFS_134samples"
cd $wrkDir/$subDir
inFile="allsites.25samples.noORF.noMiss.500.SNVs.snpden.chr"
vcfNoMissingOnlySNVs="vcfNoORFsSNVs.vcf"
nbOfChroms=16

#chrNb=$1

# output
oneSNPVCF="allChrom.1SNP.25samples.noORF.noMiss.SNVs.demoHist"

# check if the input vcf is indexed
gzVCFNoMisOnlySNVs="allsites.25samples.noORF.noMiss.SNVs.demoHist.vcf.gz"
if [ ! -f "$gzVCFNoMisOnlySNVs.tbi" ];
then
    bcftools index -t $gzVCFNoMisOnlySNVs
fi

for ((chrNb = 1 ; chrNb <= $nbOfChroms ; chrNb++)); 
do
    # checking if there is input vcf and chrom-specific file with 500bp windows. 
    if [ ! -f "$inFile$chrNb" ] | [ ! -f "$vcfNoMissingOnlySNVs" ]; 
    then
    echo "ERROR: no inputfiles - $inFile$chrNb or $vcfNoMissingOnlySNVs - found. "
    exit 0 
    fi

    #head -20 $inFile$chrNb

    # the command line below starts with the first window and keeps those windows with SNPs ($4>1)
    # and 1000 apart from each other. 
    beg=0; awk -v b=$beg '{if($4>=1 && $2>=b) {print $0; b = $3+1000}}' allsites.25samples.noORF.noMiss.500.SNVs.snpden.chr$chrNb \
    > allsites.25samples.noORF.noMiss.500.SNVs.indepWindows.chr$chrNb

    awk 'BEGIN{srand();} 
    FNR == NR { chr[FNR] = $1; 
    start[FNR] = $2; 
    end[FNR] = $3; 
    snpCounts[FNR]=$4; 
    randSNP[FNR]=int(rand()*$4)+1; 
    intCounts[FNR]=0
    count++; next } 
    FNR != NR {
        #print $0
        for (i=1; i<=count; i++)

            if($1==chr[i] && $2>=start[i] && $2<end[i]) {
                intCounts[i]++;
                if(intCounts[i]==randSNP[i]) {
                    print chr[i],start[i],end[i],snpCounts[i],randSNP[i],$1,$2
                    break;
                }
                break;
            };       
    }
    ' allsites.25samples.noORF.noMiss.500.SNVs.indepWindows.chr$chrNb $vcfNoMissingOnlySNVs > chromosome$chrNb.1SNPperWindow.tmp
done

fCOUNT=0
for ((chrNb = 1 ; chrNb <= $nbOfChroms ; chrNb++));  
do
    if [ -f chromosome$chrNb.1SNPperWindow.tmp ];
    then
        fCOUNT=$((fCOUNT + 1))
        cat chromosome$chrNb.1SNPperWindow.tmp >> allChromosomes.1SNPperWindow.txt
    fi
done
if [ $fCOUNT -lt 16 ];
then
    echo "ERROR: not all chromosome*.1SNPperWindow.tmp were generated"
    exit 0;
else
    cat allChromosomes.1SNPperWindow.txt | awk '{print $6,$7}' | sed 's/ /\t/g' > snpsToKeep.1SNPperWindow.tmp
    bcftools view -R snpsToKeep.1SNPperWindow.tmp $gzVCFNoMisOnlySNVs -Oz -o $oneSNPVCF.vcf.gz
    #rm -f *.tmp allsites.25samples.noORF.noMiss.500.SNVs.indepWindows.chr*
    vcftools --gzvcf $oneSNPVCF.vcf.gz --geno-r2 --ld-window-bp 10000 --out rsquared.1SNPwindow
    bcftools query -f '[%GT\t]\n' $oneSNPVCF.vcf.gz > $oneSNPVCF.GT
    bcftools query -f '%CHROM\t%POS\n' $oneSNPVCF.vcf.gz > $oneSNPVCF.chrpos
    #8.2 - Re-coding GT table 
    cat $oneSNPVCF.GT | sed 's/|/\//g' | sed 's/\.\/\./-1/g' | sed 's/0\/0/0/g' | sed 's/0\/1/1/g' | sed 's/1\/0/1/g' | sed 's/1\/1/2/g' > $oneSNPVCF.coded.GT

 
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
