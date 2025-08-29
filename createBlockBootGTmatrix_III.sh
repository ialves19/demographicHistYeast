#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 8                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o bootRd.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e bootRd.%N.%j.err.log      # File to which STDERR will be written

display_help () {

    echo ""
    echo "# HELP: -h"
    echo ""
    echo "This script performs the STEP THREE of the blockbootstrap approach "
    echo "used for demographic inference - parameter estimation."
    echo "It generates a GT matrix with a SNP copied as many times as"
    cho "the number of times a genomic window has been drawn in the 1st step."
    echo ""
    echo "It requires a file containing genomic regions of an arbitrary size."
    echo "The file must have the following structure: "
    echo "<chromosome> <startCoordinate> <endCoordinate> <nb of sites surveyed (i.e, passed QC)> <nb of times sampled>"
    echo "It requires the ORIGINAL GT MATRIX with the corresponding .chrPos file."
    echo ""
    echo "Such file must be previously generated using the blockBootComputation_I.sh"
    echo "(1) The script will go through the ORIGINAL GT matrix"
    echo "(2) It will identify within each window it falls"
    echo "(3) It checks how many times such window has been bootstraped"
    echo "(4) and copies a SNP this amount of times"
    echo ""
    echo "The script should be run as: "
    echo "sbatch createBlockBootGTmatrix_III.sh -i <inputDir> -g <ORIGINAL GT matrix - NO EXT> -p <prefix of boot vcf-style GT> -b <boot nb.>"
    echo "" 

}

# creates a matrix where each SNP from the vcfLikeMatrix 
# is repeated as many times as the number of times a windows
# was picked up in the blockbootstrap
creatingGTmaBlockboot () {

    local oFile="bootstrap_DS_${bootNb}/$fbootPrefix.boot$bootNb.vcf.GT"

    awk 'FNR == NR { chr[FNR] = $1; start[FNR] = $2; end[FNR] = $3; sitesSurveyed[FNR]=$4; counts[FNR]=$5; count++; next } 
    FNR != NR {
        #print $0
        for (i=1; i<=count; i++)

            if($1==chr[i] && $2>start[i] && $2<=end[i]) {
                if(counts[i]>0) {
                    inc = 1;
                    while(inc <= counts[i]) {
                        print $0;
                        inc++;
                    };
                };
                break;
            };
    }' $bootFile $vcfLikeMatrix > $oFile
}

# verifying blockboostraps 
verifyingBlockBootstraps () {

    local tmpGTm="bootstrap_DS_${bootNb}/tmp.${bootNb}.GT"
    sort bootstrap_DS_${bootNb}/$fbootPrefix.boot$bootNb.vcf.GT | uniq -c | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' > $tmpGTm
    awk 'BEGIN { errorSum=0 } FNR == NR { chr[FNR] = $1; start[FNR] = $2; end[FNR] = $3; sitesSurveyed[FNR]=$4; counts[FNR]=$5; count++; next } 
    FNR != NR {
        #print $0
        for (i=1; i<=count; i++)

            if($2==chr[i] && $3>start[i] && $3<=end[i]) {
                
                #print chr[i],start[i],end[i],counts[i],$1,$2,$3,$4,$5,$6;
                if(counts[i]!=$1) {
                    print chr[i],start[i],end[i],counts[i],$1,$2,$3,$4,$5,$6;
                    errorSum++;
                }
                break;
            };          
    }
    END {
        print errorSum;
    }' $bootFile $tmpGTm
    echo "Verification DONE";
    echo "Cleaning up";

    rm -f $tmpGTm
}
######################################
##
##           END of FUNCTIONS
##
######################################

######################################
##
##           MAIN 
##
######################################

start=`date +%s`


if [ "$1" == "-h" ]; 
then
    display_help;
    echo ""
else

    while getopts "i:g:p:b:" option
    do 
        case "${option}"
            in
            
            i)inputDir=${OPTARG};;
            g)inGTmatrix_AA=${OPTARG};;
            p)fbootPrefix=${OPTARG};;
            b)bootNb=${OPTARG};;
            
        esac
    done

    # Debug
    # inputDir="/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples"
    # inGTmatrix_AA="samples25.coded.AA.all"
    # fbootPrefix="allsites.25samples.noORF.noMiss.500.siteden"
    # bootNb=1

    echo "Working on dir: $inputDir";
    echo "Working on GT matrix: $inGTmatrix_AA.GT";
    echo "and the corresponding SNP pos: $inGTmatrix_AA.chrPos";
    echo "Bootstrap: $bootNb";
    echo "";

    #----------------------
    # Checking if dir/files exist
    #----------------------
    FCOUNT=0
    if [ ! -f "$inputDir/$inGTmatrix_AA.GT" ] | [ ! -f "$inputDir/$inGTmatrix_AA.chrPos" ]; 
    then
        echo "ERROR: no GT nor SNP matrix. "
        exit 1
    else 
        FCOUNT=$((FCOUNT+1))
    fi

        #----------------------
    # Blockbootstrap computation
    # main 
    #----------------------
    if [ $FCOUNT -gt 0 ]
    then
        cd $inputDir
        bootFile="bootstrap_DS_${bootNb}/$fbootPrefix.boot$bootNb.fullDS"
        vcfLikeMatrix="$inGTmatrix_AA.vcf.GT"
        paste -d$'\t' $inGTmatrix_AA.chrPos $inGTmatrix_AA.GT > $vcfLikeMatrix
        
        if [ -f $bootFile ]
        then
            creatingGTmaBlockboot;
            blockBootGTfile="bootstrap_DS_${bootNb}/$fbootPrefix.boot$bootNb.vcf.GT"
            #verifyingBlockBootstraps;

        else 
            echo "ERROR: no bootstrap file found. "
            exit 1
        fi
    fi
fi
 
echo "DONE!";
echo "";

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