while read -r chrnb pos nbsnps densnps; do
   upBound=$(( $pos + 500 )) 
   echo "Processing window: $pos $upBound";
   nbBPs=`less $outputVCForf | grep $chrnb$'\t' | awk -v beg=$pos -v end=$upBound '{if ($2 >= beg && $2 < end) print $0; else if ($2 > end) exit 0; }' | wc -l`
   echo "$chrnb $pos    $upBound   $nbBPs" >> nbBasePairPerWindowCOCO.txt
done < noheader.allsites.25samples.0.5kb.snpden


while read -r chrnb pos nbsnps densnps;
do
    upBound=$(( $pos + 500 ))
    echo "Processing window: $chrnb $pos $upBound";
    #nbBPs=`less $outputVCForf | grep $chrnb$'\t' | awk -v beg=$pos -v end=$upBound '{if ($2 >= beg && $2 < end) print $0; else if ($2 > end) exit 0; }' | wc -l`
    less $outputVCForf |  grep $chrnb$'\t' | while read line;
    do 
        chrPos=`cut -d$'\t' -f2 $line`; 
        if [[ $chrPos -ge $pos ]] && [[ $chrPos -lt $upBound ]] 
        then 
            echo $line; 
        elif [[ $chrPos -gt $end ]]; 
        then
            break; 
        fi; 
    done;
    #echo "$chrnb $pos    $upBound   $nbBPs" >> nbBasePairPerWindowCOCO.txt
done < noheader.allsites.25samples.0.5kb.snpden

while read -r chrnb pos nbsnps densnps;
do
    upBound=$(( $pos + 500 ))
    echo "Processing window: $chrnb $pos $upBound";
    #nbBPs=`less $outputVCForf | grep $chrnb$'\t' | awk -v beg=$pos -v end=$upBound '{if ($2 >= beg && $2 < end) print $0; else if ($2 > end) exit 0; }' | wc -l`
    less $outputVCForf |  grep $chrnb$'\t' | while read line;
    do 
        chrPos=`echo $line | cut -d$' ' -f2;`
        #chrPos=`cut -d$' ' -f2 $line`; 
        if [[ $chrPos -ge $pos ]] && [[ $chrPos -lt $upBound ]] 
        then 
            echo $line; 
        elif [[ $chrPos -gt $upBound ]]; 
        then
            break; 
        fi;
    done
done < noheader.allsites.25samples.0.5kb.snpden