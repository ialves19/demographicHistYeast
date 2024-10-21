#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o demoIII.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e demoIII.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

inDir="/shared/home/ialves/ancestralState/03-Data/allSites"
wrkDir="/shared/home/ialves/3039Sace/03-data"

orf="orf_coding.fasta"

cd $wrkDir

grep "Chr" orf_coding.fasta | cut -d$',' -f2 | sed 's/from //g' | \
sed 's/-/ /g' | sed 's/^ //g' | sed 's/ /\t/g' | sed 's/Chr/chromosome/g' | grep -v "Mito" > ORF.Sace.coordinates.tmp

paste <(cut -d$'\t' -f1 ORF.Sace.coordinates.tmp) <(cut -d$'\t' -f2 ORF.Sace.coordinates.tmp) | \
    sed 's/\t//g' > chrID.tmp
paste <(cut -d$'\t' -f3-4 ORF.Sace.coordinates.tmp) <(cat chrID.tmp) > ORF.Sace.coordinates.roman
rm -f ORF.Sace.coordinates.tmp
chrNbRoman=( "I" "II" "III" "IV" "V" "VI" "VII" "VIII" "IX" "X" "XI" "XII" "XIII" "XIV" "XV" "XVI" ) 
COUNT=1
for chrID in ${chrNbRoman[@]}; 
    do 
        echo "Chromosome $chrID has: "
        grep chromosome$chrID$ ORF.Sace.coordinates.roman | wc -l
        
        sed -i "s/chromosome$chrID$/chromosome$COUNT/g" ORF.Sace.coordinates.roman 
        cat ORF.Sace.coordinates.roman | grep "chromosome$COUNT$" > ORF.Sace.coordinates.chromosome$COUNT.tmp
        # because some regions are reverse complement one need to reverse hte
        # beginning and the end of the region. 
        cat ORF.Sace.coordinates.chromosome$COUNT.tmp | awk '{OFS="\t"} {if ( $1 > $2 ) {print $3,$2,$1} else {print $3,$1,$2} }' > ORF.Sace.coordinates.chromosome$COUNT

        ((COUNT++));
done

rm -f ORF.*.tmp
COUNT=1; 
for file in ORF.Sace.coordinates.chromosome{1..16}; 
do 
    if [[ $COUNT -eq 1 ]]; 
    then 
        cat $file > ORF.Sace.coordinates.allChroms; 
    else 
        cat $file >> ORF.Sace.coordinates.allChroms;
    fi; 
    ((COUNT++)); 
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
