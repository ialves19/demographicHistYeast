#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 12                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=500              # Memory per cpu
#SBATCH -o demoIV.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e demoIV.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

wrkDir="/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/input/TwoDomest_EurOld_dSFS"
cd $wrkDir

fscFolder="fsc27_linux64"

fsc="fsc27093"
cp ../../${fscFolder}/$fsc .

jobcount=0
msgs=conOutputs

#-------- Number of different runs per data set ------
numRuns=1
runBase=1 
#-----------------------------

mkdir $msgs 2>/dev/null

#-------- Default run values ------
numSims=100000                    #-n command line option
numCycles=40                      #-L command line option
minValidSFSEntry=1                #-C command line option

#-------- Ascertainment ------
#withAscertainment=0
#ascPop=0                          #-a command line option
#ascSize=2                         #-A command line option
#-----------------------------
useMonoSites="" 
#useMonoSites="-0"                #-0 command line option
#----------multiSF------------
#multiSFS=""
multiSFS="--multiSFS"            #--multiSFS command line option
#-----------------------------

#-------- Generic Name ------
genericName="OneDomest_3_1"
tplGenericName="OneDomest_3_1"
estGenericName="OneDomest_3_1"
echo "Running fastsimcoal for model: $genericName" 
#-----------------------------

./$fsc -t ${tplGenericName}.tpl -e ${estGenericName}.est -n ${numSims} -M -L ${numCycles} -d -C 1 -c0 -B12 $multiSFS -q



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
