#!/bin/bash

#SBATCH -p public                 # Partition to submit to (fast / slow)
#SBATCH -N 1                      # Nb of nodes
#SBATCH -c 1                      # Nb of CPUs per task 
#SBATCH -n 1                      # Nb tasks per CPU
#SBATCH --job-name="collectFsc"              # Memory per cpu
#SBATCH -o testScal.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e testScal.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

wrkDir="/home2020/home/gmgm/ialves/demoHist_yeast3039/04-analysis/fastsimcoal/input/SCALLING"
cd $wrkDir

fscFolder="fsc27_linux64"

fsc="fsc27093"
#cp ../../${fscFolder}/$fsc .

jobcount=0
msgs=conOutputs

#-------- Number of different runs per data set ------
numRuns=10
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
genericName="OneDomest"
cd ${genericName}_dSFS/ 

echo "Collecting fastsimcoal runs for model: $genericName" 
#-----------------------------
for (( runsDone=$runBase; runsDone<=$numRuns; runsDone++ ))
do
        runDir="${genericName}_$runsDone"
        echo "Extracting parameters from run $runsDone"
        #Processing best likelihood files
        bestlhoodFile=${runDir}.bestlhoods
        #Extract second line
        if [ $runsDone -eq 1 ];
                then
                        header=$(sed '1!d'  ${runDir}/${runDir}/$bestlhoodFile)
                        echo -e "Run\t$header" > ${genericName}_bestLhoodsParams.txt
        fi
        wantedParameters=$(sed '2!d'  ${runDir}/${runDir}/$bestlhoodFile)
        echo -e "$runDir\t$wantedParameters" >> ${genericName}_bestLhoodsParams.txt
done


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