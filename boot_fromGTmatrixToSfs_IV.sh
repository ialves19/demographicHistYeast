#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o bootIV.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e bootIV.%N.%j.err.log      # File to which STDERR will be written

display_help () {

    echo ""
    echo "# HELP: -h"
    echo ""
    echo "This script performs the STEP FOUR of the blockbootstrap approach "
    echo "used for demographic inference - parameter estimation."
    echo "It generates a 1D-SFS, 2D-SFS and multi-SFS from the GT matrix generated under"
    echo "each bootstrap as in step III. "
    echo ""
    echo "It requires as input files:" 
    echo "-s <A file with samples ordered as in the genotype matrix (-g); ex: samples25.txt>"
    echo "-d. <A list of samples with the respective clade; ex: list_of_strains_per_Clade.txt;"
    echo "-z. df_zygosity_allHet.txt; file which defines whether samples should be treated as homozygous or heterozygous."
    echo "-g.  genotype matrix generated in previous step; ex: allsites.25samples.noORF.noMiss.500.siteden.boot1.vcf.GT;"
    echo "-b T/F; information to pass to the script computingOBS_1D_SFS_6.sh which will lead a change in the output folder and re-direct the output to the bootDir. "
    echo "-n 1..100; bootstrap number to be taken into account in the bootDir and output files. "
    echo "-o <output directory where the folder containing the 2D/multi-SFS will be generated (NEEDED FOR THE D/multi-SFS)>; ex: /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/fiveClades"
    echo "-m <file containing the name of the model that will give rise to the folder. It also contains the order of the samples in the model and the zygosity. (NEEDED FOR THE 2D/multi-SFS); ex: model_TwoDomest_EurOld_Taiw" 
    echo ""
    echo "The script should be run as: "
    echo "sbatch boot_fromGTmatrixToSfs_IV.sh -d /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples -s samples25.txt -f list_of_strains_per_Clade.txt -z df_zygosity_allHet.txt -g allsites.25samples.noORF.noMiss.500.siteden.boot1.vcf.GT -b T -n 1 -o /shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/fiveClades -m model_TwoDomest_EurOld_Taiw_stru.txt"
    echo "" 
    echo "by Isabel Alves - June 2025" 

}

if [ "$1" == "-h" ]; 
then
    display_help;
    echo ""
else
    start=`date +%s`
    conda activate rForDemoInf

    while getopts "d:s:f:z:g:b:n:o:m:" option
    do 
        case "${option}"
            in
            
            d)wrkDir=${OPTARG};;
            s)samples=${OPTARG};;
            f)dfstrains=${OPTARG};;
            z)zig=${OPTARG};;
            g)gtmatrix=${OPTARG};;
            b)boot=${OPTARG};;
            n)bootNb=${OPTARG};;
            o)outDir=${OPTARG};;
            m)model=${OPTARG};;
            
        esac
    done

    ## debugging
    #wrkDir="/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/SFS_134samples"
    #outDir="/shared/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/fiveClades"
    #samples="samples25.txt"
    #dfstrains="list_of_strains_per_Clade.txt"
    #zig="df_zygosity_allHet.txt"
    #gtmatrix="allsites.25samples.noORF.noMiss.500.siteden.boot1.vcf.GT"
    #model="model_TwoDomest_EurOld_Taiw_stru.txt"
    #boot=T
    #bootNb=1

    echo "Converting GT matrices to 1DSFS. Working dir: $wrkDir"
    echo "Generating GT matrices to 2D/multi-SFS. Model: $model"
    echo "Remaining files can be found : $outDir"
    echo ""
    echo "Input files and variables:"
    echo "Matrix of Positions: $gtmatrix";
    echo "List Samples: $samples";
    echo "Strain df: $dfstrains";
    echo "Zygosity df: $zig";
    echo "Analysing bootstraps? $boot";
    echo "Bootstrap nb: $bootNb";
    echo ""

    #----------------------
    # Checking if dir/files exist
    #----------------------
    FCOUNT=0
    if [ ! -f "$wrkDir/bootstrap_DS_${bootNb}/$gtmatrix" ]; 
    then
        echo "ERROR: no GT matrix. "
        exit 1
    elif [ ! -f "$wrkDir/$samples" ];
    then
        echo "ERROR: no file with sample order. "
        exit 1        
    elif [ ! -f "$wrkDir/$dfstrains" ];
    then 
        echo "ERROR: no file with clade/strain correspondance. "
        exit 1
    elif [ ! -f "$wrkDir/$zig" ];
    then 
        echo "ERROR: no file with zygosity for 1D-sfs computation "
        exit 1
    elif [ ! -f "$outDir/$model" ];
    then 
        echo "ERROR: no file with model specification for 2D-sfs computation "
        exit 1
    else
    FCOUNT=$((FCOUNT+1));
    fi

    if [ $FCOUNT -gt 0 ]
    then
        cd $wrkDir

        # creating intermediary file 
        cut -d$'\t' -f8- bootstrap_DS_$bootNb/$gtmatrix > bootstrap_DS_$bootNb/tmp.GT

        cd ~/demoHist_yeast3039/02-scripts
        chmod +x computingOBS_1D_SFS_6.R
        chmod +x computingOBS_2D_multi_SFS_7.R
        chmod +x correctingMONOsitesBoot.R

        echo "Computing 1D-SFS...";
        sbatch computingOBS_1D_SFS_6.sh -d $wrkDir -s $samples -f $dfstrains -z $zig -g bootstrap_DS_$bootNb/tmp.GT -b $boot
        echo "";
        echo "Computing 2D-SFS and multi-SFS...";
        sbatch computingOBS_2D_multi_SFS_7.sh -i $wrkDir -o $outDir -s $samples -d $dfstrains -m $model -g bootstrap_DS_$bootNb/tmp.GT -b $boot
        modelTag=`echo $model | sed 's/model_\(.*\).txt/\1/'`
        while true;
        do
            if [ -f "$outDir/${modelTag}_${boot}_dSFS/${modelTag}_${boot}_DSFS.tmp" ];
            then 
                echo "sfs.tmp found. Correcting number of mono sites."
                Rscript correctingMONOsitesBoot.R $wrkDir $outDir $modelTag $boot
                break;
            fi
        done
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

fi


