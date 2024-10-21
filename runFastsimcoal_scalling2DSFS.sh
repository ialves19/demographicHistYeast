#!/bin/bash

#SBATCH -p public                 # Partition to submit to (fast / slow)
#SBATCH -N 1			  # Nb of nodes
#SBATCH -c 1 			  # Nb of CPUs per task	
#SBATCH -n 1			  # Nb tasks per CPU
#SBATCH --job-name="fscoal"              # Memory per cpu
#SBATCH -o testScal.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e testScal.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

wrkDir="/home2020/home/gmgm/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/input/2D-SFS_SCALLING"
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
multiSFS=""
#multiSFS="--multiSFS"            #--multiSFS command line option
#-----------------------------

#-------- Generic Name ------
genericName="OneDomestTaiw"
cd ${genericName}_dSFS/
mkdir ${genericName}_${1}
cd ${genericName}_${1}/
cp ../../../../${fscFolder}/$fsc .
cp ../${genericName}_jointDAFpop* .
for file in ${genericName}_jointDAFpop*; 
do 
	sufix=`echo $file | sed 's/.*\(_joint.*\)/\1/'`; 
	cp $file ${genericName}_${1}${sufix};
done
#cp ../${genericName}_DSFS.obs ${genericName}_${1}_DSFS.obs 
cp ../${genericName}.tpl ${genericName}_${1}.tpl
cp ../${genericName}.est ${genericName}_${1}.est
tplGenericName="${genericName}_${1}"
estGenericName="${genericName}_${1}"

echo "Running fastsimcoal for model: $genericName job number ${1}" 
#-----------------------------

./$fsc -t ${tplGenericName}.tpl -e ${estGenericName}.est -n ${numSims} -M -L ${numCycles} -d -C1 -c1 -B1 $multiSFS -q


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
