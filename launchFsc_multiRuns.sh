#!/bin/bash
#SBATCH -p gpu                 # Partition to submit to (fast / slow)
#SBATCH -n 10                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=500              # Memory per cpu
#SBATCH -o test.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e test.%N.%j.err.log      # File to which STDERR will be written

numRuns=50
#-------- Working directory ------
wrkDir="/home/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/input"
cd $wrkDir
#-------- fsc directory ------
fscFolder="fsc27_linux64"
#-------- fsc exec ------
fsc="fsc27093"
#-------- Generic Name ------
genericName="OneDomest_dSFS"
cd $genericName
#-------- Default run values ------
numSims=100000                   #-n command line option
numCycles=40                      #-L command line option
minValidSFSEntry="-C1"                #-C command line option
#--------- More fsc settings --------
# !!! comment out if needed !!!
useMonoSites="" 
#useMonoSites="-0"                #-0 command line option
#multiSFS=""
multiSFS="--multiSFS"            #--multiSFS command line option
#-----------------------------

# get through all your fastq files
for (( i=2; i<${numRuns}; i++ ))
do

  mkdir ${genericName}_${i}
  cp ../../${fscFolder}/$fsc ${genericName}_${i}/
  cp ${genericName}.est ${genericName}_${i}/${genericName}_${i}.est
  cp ${genericName}.tpl ${genericName}_${i}/${genericName}_${i}.tpl
  cp ${genericName}_DSFS.obs ${genericName}_${i}/${genericName}_${i}_DSFS.obs
  cd ${genericName}_${i}/
  echo "Running model: ${genericName}, independent run ${i}"
  srun --exclusive -N1 -n1 ./$fsc -t ${genericName}_${i}.tpl -e ${genericName}_${i}.est -n ${numSims} -M -L ${numCycles} -d ${minValidSFSEntry} -c ${SLURM_CPUS_PER_TASK} -B12 $multiSFS -q &
  cd ..

done
wait
