#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -c 1				# Nb of tasks
#SBATCH -c 8    			# Nb of CPUs per task
#SBATCH --job-name="fscoal"              # Memory per cpu
#SBATCH -o testScratch.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e testScratch.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

wrkDir="/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/input"
cd $wrkDir

fscFolder="fsc27_linux64"

fsc="fsc27093"
#cp ../../${fscFolder}/$fsc .

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
genericName="TwoDomest_EurOld"
cd ${genericName}_dSFS/
mkdir ${genericName}_${1}
cd ${genericName}_${1}/
cp ../../../${fscFolder}/$fsc .
cp ../${genericName}_DSFS.obs ${genericName}_${1}_DSFS.obs 
cp ../${genericName}.tpl ${genericName}_${1}.tpl
cp ../${genericName}.est ${genericName}_${1}.est
tplGenericName="${genericName}_${1}"
estGenericName="${genericName}_${1}"

echo "Running fastsimcoal for model: $genericName job number ${1}" 
#-----------------------------

./$fsc -t ${tplGenericName}.tpl -e ${estGenericName}.est -n ${numSims} -M -L ${numCycles} -d -C1 -c8 -B8 $multiSFS -q


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
