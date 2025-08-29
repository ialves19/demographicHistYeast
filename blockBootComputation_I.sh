#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 8                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o bootst.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e bootst.%N.%j.err.log      # File to which STDERR will be written

display_help () {

    echo ""
    echo "# HELP: -h"
    echo ""
    echo "This script performs the STEP ONE of the blockbootstrap approach used for demographic inference - parameter estimation."
    echo "It randomly samples genomic windows and creates ONE bootstrap from a original genomic dataset. "
    echo ""
    echo "It requires a file containing genomic regions of an arbitrary size."
    echo "The file must have the following structure: "
    echo "<chromosome> <startCoordinate> <endCoordinate> <nb of sites surveyed (i.e, passed QC)>"
    echo ""
    echo "Such file must be previously generated (see gettingSNPdensity.bash)"
    echo "(1) The script will randomly select windows"
    echo "(2) It will count how many sites have been surveyed after resampling"
    echo "(3) It will compare whether the bootstrap is +/- 1% of the full dataset in terms of total nb of sites surveyed."
    echo "(4) It will compute a block-bootstrap following the 1% condition"
    echo "(5) It will generate a new file with the exact same structure as the input file but "
    echo "with an extra column (5th column) indicating how many times a window has been sampled."
    echo ""
    echo "The script should be run as: "
    echo "sbatch blockBootComputation_I.sh -i <inputDir> -d <input Site Den> -b <bootNb>"
    echo "" 

}

# creating bootstrap folder 
creatingSubDirBoot () {

    if [ ! -d "$inputDir/bootstrap_DS_${bootNb}" ];
    then 
        echo "Creating bootstrap folder nb. $bootNb"
        mkdir bootstrap_DS_${bootNb}
        #mkdir bootstrap_DS_${bootNb}/tmp_${bootNb}_${chromNb}
    else 
        echo "Bootstrap folder nb. $bootNb exists"
        
        # if [ -d "bootstrap_DS_${bootNb}/tmp_${bootNb}_${chromNb}" ];
        # then 
        #     echo "Deleting previous tmp_${bootNb}_${chromNb}"
        #     rm -rf bootstrap_DS_${bootNb}/tmp_${bootNb}_${chromNb}
        #     echo "Creating a new tmp_${bootNb}_${chromNb}"
        #     mkdir bootstrap_DS_${bootNb}/tmp_${bootNb}_${chromNb}
        # else 
        #     mkdir bootstrap_DS_${bootNb}/tmp_${bootNb}_${chromNb}
        # fi
    fi
    bootFullPath="bootstrap_DS_${bootNb}"
    echo $bootFullPath;
}
# creating tmp_1_1 folder where 1_1 means 1 bootstrap 1 chromosome
creatingTmpDirBoot () {
    
    local chrNb=$1
    if [ ! -d "$bootFullPath/tmp_${bootNb}_${chrNb}" ];
    then 
        echo "Creating tmp folder nb. $bootNb, chrom. ${chrNb}"
        mkdir $bootFullPath/tmp_${bootNb}_${chrNb}
    else 
        echo "Temporary folder: $bootFullPath/tmp_${bootNb}_${chrNb} exists"
        echo "Deleting previous tmp_${bootNb}_${chrNb}"
        rm -rf $bootFullPath/tmp_${bootNb}_${chrNb}
        mkdir $bootFullPath/tmp_${bootNb}_${chrNb}

    fi

    tmpDirBoot="$bootFullPath/tmp_${bootNb}_${chrNb}";
}
# count numbers in the forth column of a matrix
countNbColFour () {

    local fName=$1
    local sumCol4=`cat $fName | awk '{sum+=$4}END{print sum}'`
    echo "$sumCol4"

}

# get nb of Windows and nb of surveyed sites across windows
countWindowsAndSurveyedSites () {

    local fName=$1
    local tmpFolder=$2
    local chrNb=$3
    # selecting chromosome 
    cat $fName | grep chromosome$chrNb$' ' > $tmpFolder/$fName.chr${chrNb}
            
    # getting the number of windows for the chromosome
    nbWindows=`cat $tmpFolder/$fName.chr${chrNb} | wc -l`
    # number of surveyed sites
    # surveyed sites per window are in column nb 4
    nbOfSeqSites="$(countNbColFour "$tmpFolder/$fName.chr${chrNb}")"
    #echo $nbOfSeqSites

}

# re-sample with replacement from a file with X lines (i.e windows or block)
randomlySelectWindows () {

    local fName=$1
    
    shuf -r -n $nbWindows $fName > $fName.boot${bootNb}
    local bootTmpFName="$fName.boot${bootNb}"
    nbOfSeqSitesBoot="$(countNbColFour $bootTmpFName)"

}

# re-sample with replacement from a file with X lines (i.e windows or block)
# until the total nb of surveyed sites is within 1% of that of the ORIGINAL DATASET
blockBootstrap () {

    local x1=$1 # nb of Surveyed Sites
    # if the amount of surveyed sites (mono+polymorphic) in the bootstrap is >  than 1% the resampling is done 
    while true;
    do
        # randomly drawing windows
        randomlySelectWindows $chrFileName 
        # checking how many sites of the randomly selected are present in the original dataset. 
        ratioBoot=`awk -v "x=$x1" -v "y=$nbOfSeqSitesBoot" 'BEGIN {z=1-(x/y); if(z < 0.0) {print -z} else {print z}}'`

        # bash does not handle non-integers.
        # we need to use following trick
        comparison=$(echo "$ratioBoot > $theshold" | bc)
        if [[ $comparison -eq 1 ]];
        then
            #echo $ratioBoot
            rm -f $chrFileName.boot${bootNb}
        else
            break;
        fi
    done
    # END of the resampling 
} 
# END of the block bootstrap function 

# re-sample with replacement from a file with X lines (i.e windows or block)
# until the total nb of surveyed sites is within 1% of that of the ORIGINAL DATASET
annotatingWindowFileWithBootstraps () {

local origFile=$1 
local bootFileSort=$2
# below we open the two files $inputSiteDen.chr${chromNb} 
# $inputSiteDen.chr${chromNb}.boot${bootNb}.sorted
# with the first one chrX.siteden in memory we go over all the positions within the resampled file
# we count how many times each window was resampled in a given bootstrap and add +1
# this gives the nb of times a window of the original dataset was resampled. 
awk 'FNR == NR { chr[FNR] = $1; start[FNR] = $2; end[FNR] = $3; sitesSurveyed[FNR]=$4; counts[FNR]=0; count++; next } 
FNR != NR {
    #print $0
    for (i=1; i<=count; i++)

        if($1==chr[i] && $2==start[i] && $3==end[i]) {
            counts[i]++;
            #print chr[i],start[i],end[i],counts[i];
            break;
        };          
}
END {
    for (i=1; i<=count; i++)
        print chr[i],start[i],end[i],sitesSurveyed[i],counts[i] 
}' $origFile $bootFileSort > $bootFullPath/$inputSiteDen.boot${bootNb}.chr${chromNb}

}
######################################
##
##
##      END OF FUNCTIONS
##
##
######################################

######################################
##
## This script creates X bootstraps from the original dataset for 25 or
## 27 samples using the Taiwan Wild 1 or the Wild Chinese as a wild group. 
## 
## This script takes the files containing the information about 
## the nb of sites sequenced within XXbp windows. 
## Files: allsites.25samples.noORF.noMiss.500.siteden 
## 
## These files can be found across the multiple folders: 
## /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_*samples
## within pangloss. 
## The file *.siteden was generated with the script: gettingSNPs_density. 
## input files: 
## 1. allsites.25samples.noORF.noMiss.500.siteden 
## 2. ORF.Sace.coordinates.allChroms - file should be generated with preparingORFregions.sh
## 
## output files: 
## 1. allsites.25samples.noORF.noMiss.SNVs.demoHist.vcf.gz/vcfNoORFsSNVs.vcf
## 2. allsites.25samples.noORF.noMiss.500.SNVs.snpden.chr{1..16}
## the (1) is the clean vcf : no missing data and only SNVs
## the (2) is the list of 500bp (TO BE SPECIFIED) regions with the number of SNVs per window
## (2) is the input of XXXXXXX.sh which randomly selects 1SNP/window
## across windows 1000bp apart (TO BE SPECIFIED)
## 3. The files (2) are merged into a single file *boot*.fullDS
##
##
## by Isabel Alves - April 2025
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

    while getopts "i:d:b:" option
    do 
        case "${option}"
            in
            
            i)inputDir=${OPTARG};;
            d)inputSiteDen=${OPTARG};;
            b)bootNb=${OPTARG};;
            
        esac
    done

    # Debug
    #inputDir="/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples"
    #inputSiteDen="allsites.25samples.noORF.noMiss.500.siteden"
    #chromNb=1
    #bootNb=1

    echo "Working on dir: $inputDir";
    echo "Working on file: $inputSiteDen";
    echo "Bootstrap: $bootNb";
    echo "";

    #----------------------
    # Setting variables
    #----------------------
    # the bootstrap should contain around +- 1% of the total number of surveyed sites 
    theshold=0.01
    echo "Bootstrap datasets must contain +/- : $theshold of the sites in the original dataset";
    echo "";
    nbChroms=16

    #----------------------
    # Checking if dir/files exist
    #----------------------
    FCOUNT=0
    if [ ! -f "$inputDir/$inputSiteDen" ] | [ ! -d "$inputDir" ]; 
    then
    echo "ERROR: no inputfile or no input dir. "
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
        # creating folder structure - boot folder
        creatingSubDirBoot;
        for ((chromNb=1; chromNb<=nbChroms; chromNb++))
        do
            echo "**Chromosome: $chromNb **";
            echo "";
            if [ ! -z "$bootFullPath" ]
            then 
                # creating folder structure - tmp folder
                creatingTmpDirBoot $chromNb;
                echo "Creating boot nb: $bootNb for chromosome nb: $chromNb within: $tmpDirBoot"
                # getting nb of windows and total nb of surveyed sites ORIGINAL DATASET
                countWindowsAndSurveyedSites $inputSiteDen $tmpDirBoot $chromNb
                chrFileName="$tmpDirBoot/$inputSiteDen.chr${chromNb}"
                echo "Chromosome $chromNb has $nbWindows windows and $nbOfSeqSites surveyed sites."
                # creating a blockbootstrap containing +/- 1% of the surveyed sites in ORIGINAL DATASET
                blockBootstrap $nbOfSeqSites
                # randomly selecting windows with replacement 
                echo "Total number of sites within bootstrap ${bootNb} is: ${nbOfSeqSitesBoot}"
                # sorting 
                cat $chrFileName.boot${bootNb} | sort -n -k2 > $chrFileName.boot${bootNb}.sorted
                rm -f $chrFileName.boot${bootNb}
                # adding a 5th column to the ORIGINAL DATASET with the nb of times a windows was resampled
                annotatingWindowFileWithBootstraps $chrFileName $chrFileName.boot${bootNb}.sorted
                
            else
                echo "ERROR: folder arborecence messed up"
                exit 0 
            fi
            echo "";
            echo "Cleaning up!";
            rm -fr $tmpDirBoot 
            echo "";
        done
        totalNbOutputs=`ls $bootFullPath/$inputSiteDen.boot${bootNb}.chr* | wc -l`
        if [ $totalNbOutputs -eq $nbChroms ];
        then
            for ((chromNb=1; chromNb<=nbChroms; chromNb++));
            do 
                cat $bootFullPath/$inputSiteDen.boot${bootNb}.chr${chromNb} >> $bootFullPath/$inputSiteDen.boot${bootNb}.fullDS
            done
        else
            echo "ERROR: bootstrap files for all chromosomes are not available."
            exit 1; 
        fi
    fi
fi 

echo "Final clean up";
rm -f $bootFullPath/$inputSiteDen.boot${bootNb}.chr*

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
