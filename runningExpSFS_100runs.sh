#!/bin/bash

#SBATCH -p public                 # Partition to submit to (fast / slow)
#SBATCH -N 1			  # Nb of nodes
#SBATCH -c 1 			  # Nb of CPUs per task	
#SBATCH -n 1			  # Nb tasks per CPU
#SBATCH --job-name="maxL"              # Memory per cpu
#SBATCH -o maxL.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e maxL.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

#####################
##
##
## This script should be run:
## sbatch runningExpSFS_100runs.sh <model name> <max Lhood run> <iteration>
## example: sbatch runningExpSFS_100runs.sh OneDomest 4 1

model=$1
maxLhoodRun=$2
runNb=$3

wrkDir="/home2020/home/gmgm/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/input/SCALLING_fiveClades"
cd $wrkDir

echo "Running model: $model maxLhood run $maxLhoodRun instance $runNb"
echo "Running in folder: $wrkDir"

cd "${model}_dSFS/${model}_${maxLhoodRun}"
mkdir "${model}_${maxLhoodRun}_maxL_${runNb}"

maxLFile="${model}_${maxLhoodRun}_maxL.par"
prefixTmp=${maxLFile%.*}
cp ${model}_${maxLhoodRun}/${maxLFile} ${model}_${maxLhoodRun}_maxL_${runNb}/${prefixTmp}_${runNb}.par
cp ${model}_${maxLhoodRun}_DSFS.obs ${model}_${maxLhoodRun}_maxL_${runNb}/${prefixTmp}_${runNb}_DSFS.obs

cd ${model}_${maxLhoodRun}_maxL_${runNb}

../fsc27093 -i ${prefixTmp}_${runNb}.par -n 200000 -d --multiSFS

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
