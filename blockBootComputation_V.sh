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

start=`date +%s`
conda activate rForDemoInf