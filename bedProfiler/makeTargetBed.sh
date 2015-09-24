#!/bin/bash

docstring="DESCRIPTION \\n
Generate target bed file for which features the profile is to be created. \\n
Input features extended left and right by given amount. Extended intervals are \\n
divided in bins of given size are renamed to have the central bin as 0. \\n
\\n
Example: \\n
echo 'chr1 1000 1010' | tr ' ' '\\\\t' \ \\n
| makeTargetBed.py -i - -g genome.txt -s 3 -c \\n
\\n 
chr1    1002    1003    -3 \\n
chr1    1003    1004    -2 \\n
chr1    1004    1005    -1 \\n
chr1    1005    1006    0 \\n
chr1    1006    1007    1 \\n
chr1    1007    1008    2 \\n
chr1    1008    1009    3 \\n
"

if [[ $1 == "-h" || $1 == "--help" ]]
then
    echo -e $docstring
    exit 1
fi

while [[ $# > 0 ]]
do
key="$1"

case $key in

    -i|--inbed)
    INBED="$2"
    shift # past argument
    ;;
    -g|--genome)
    GENOME="$2"
    shift # past argument
    ;;
    -s|--slop)
    SLOP="$2"
    shift # past argument
    ;;
    -b|--binSize)
    BIN_SIZE="$2"
    shift # past argument
    ;;
    -c|--default)
    CENTER="y"
    ;;
    *)
    # unknown option
    ;;
esac
shift # past argument or value
done

if [[ $INBED == '' ]];
then 
    echo "Need input bed"
    exit 1
fi

if [[ $GENOME == '' ]];
then 
    echo "Need genome file chrom<tab><size>"
    exit 1
fi

if [[ $SLOP == '' ]];
then 
    SLOP=1000
fi

if [[ $BIN_SIZE == '' ]];
then 
    BIN_SIZE=1
fi

cat $INBED \
| python -c "
import sys
for line in sys.stdin:
    line= line.strip().split('\t')
    start= int(line[1])
    end= int(line[2])
    mid= int(round((end - start)/2))
    line[1]= start + mid
    line[2]= line[1] + 1
    print '\t'.join([line[0], str(line[1]), str(line[2])])" \
| slopBed -g $GENOME -b $SLOP \
| windowMaker -b - -i winnum -w $BIN_SIZE \
| awk -v OFS='\t' -v slop=$SLOP '{print $1, $2, $3, $4 - (slop+1)}' \
| sort -k1,1 -k2,2n -s -S 2G \
| groupBy -g 1,2,3 -c 4 -o collapse \
| python -c "
import sys
import random
for line in sys.stdin:
    line= line.strip().split('\t')
    bins= line[3].split(',')
    random.seed(line[3]) # Random but same results from same input.
    bin= random.choice(bins)
    line[3]= bin
    print '\t'.join(line)
"

exit 0

#slop=3000
#grep -P '^2\t|^1\t' ../../20150914_nucleo_5fc/ctcf_brain_embryo/fimo_out/fimo_global.bed \
#| awk '{printf "%s\t%.0f\n", $1, $2+(($3-$2)/2)}' \
#| awk '{print $0"\t"$2+1}' \
#| slopBed -g genome.tmp -b $slop \
#| windowMaker -b - -i winnum -w 1 \
#| awk -v OFS='\t' -v slop=$slop '{print $1, $2, $3, $4 - (slop+1)}' \
#| sort -k1,1 -k2,2n -k4,4R -S 2G \
#| groupBy -g 1,2,3 -c 4 -o first \
#| pigz > targets.bed.gz 