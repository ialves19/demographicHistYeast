#!/bin/bash

#SBATCH -p public                 # Partition to submit to (fast / slow)
#SBATCH -N 1			  # Nb of nodes
#SBATCH -c 1 			  # Nb of CPUs per task	
#SBATCH -n 1			  # Nb tasks per CPU
#SBATCH --job-name="fscoal"              # Memory per cpu
#SBATCH -o chTDaOfull.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e chTDaOfull.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

##################
##
## sbatch runFastsimcoal_scalling.sh modelName runNb
##
##################

wrkDir="/home2020/home/gmgm/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/input/SCALLING_fiveClades"
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
genericName=$1
runNb=$2
cd ${genericName}_dSFS/
if [ ! -d ${genericName}_${runNb}/ ]
then
	mkdir ${genericName}_${runNb}
fi

cd ${genericName}_${runNb}/
cp ../../../../${fscFolder}/$fsc .
cp ../${genericName}_DSFS.obs ${genericName}_${runNb}_DSFS.obs 
cp ../${genericName}.tpl ${genericName}_${runNb}.tpl
cp ../${genericName}.est ${genericName}_${runNb}.est
tplGenericName="${genericName}_${runNb}"
estGenericName="${genericName}_${runNb}"

echo "Running in folder: $wrkDir"
echo ""
echo "Running fastsimcoal for model: $genericName job number ${runNb}" 
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
