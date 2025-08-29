#!/bin/bash

#SBATCH -p public                 # Partition to submit to (fast / slow)
#SBATCH -N 1			  # Nb of nodes
#SBATCH -c 1 			  # Nb of CPUs per task	
#SBATCH -n 1			  # Nb tasks per CPU
#SBATCH --job-name="maxL"              # Memory per cpu
#SBATCH -o colMaxL.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e colMaxL.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

#####################
##
##
## This script should be run:
## sbatch runningExpSFS_100runs.sh <model name> <max Lhood run> <iteration>
## example: sbatch runningExpSFS_100runs.sh OneDomest 4
## it assumes we have done 100 simulations under the maxLhood model
##
## Isabel Alves 
#####################
model=$1
maxLhoodRun=$2

wrkDir="/home2020/home/gmgm/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/input/SCALLING_fiveClades"
cd $wrkDir

echo "Collecting expected lhoods according to the model: $model maxLhood run $maxLhoodRun"
echo "Running on folder: $wrkDir"

cd "${model}_dSFS/${model}_${maxLhoodRun}"

for dir in ${model}_${maxLhoodRun}_maxL_{1..100}; 
do 
if [[ "$dir" == "${model}_${maxLhoodRun}_maxL_1" ]];
then
cat $dir/$dir/$dir.lhoods > ${model}_${maxLhoodRun}_maxL.lhoods
cat $dir/$dir/${dir}_DSFS.txt > ${model}_${maxLhoodRun}_maxL_DSFS.txt
else 
sed '1d' $dir/$dir/$dir.lhoods >> ${model}_${maxLhoodRun}_maxL.lhoods
tail -n +3 $dir/$dir/${dir}_DSFS.txt >> ${model}_${maxLhoodRun}_maxL_DSFS.txt
fi
#echo $dir
done

echo "Collecting 100 Exp SFS : DONE"


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