#!/bin/bash

set -e
set -o pipefail

docstring=' Get read count from BAM file.
\n\n
This script parses the output of samtools idxstats to quickly get the total read
count from a bam file.
\n\n
Requires awk and samtools on PATH.
\n\n
USAGE \n
readCountFromBam.sh <aln.bam> <aln2.bam> ... \n

\n\n
OUTPUT \n
Tab separated table with columns:\n
1. Input file name \n
2. Count mapped (sum 3rd col) \n
3. Count mate unmapped (sum 4th col) \n
4. Count unmapped (\* line) \n
5. Ref size (sum 2nd col) \n
'

bams="$*"

if [[ ${bams} = "" || ${bams} = "-h" || ${bams} == "--help" ]]
then
    echo -e $docstring
    exit 1
fi

# ----------------------------------------------------------

for bam in ${bams}
do
if [[ ! -f ${bam}.bai && ! -f ${bam%.bam}.bai ]]
then
    echo "Index file ${bam}.bai or ${bam%.bam}.bai not found"
    exit 1
fi

s=$(type -P samtools)
if [[ ${s} = "" ]]
then
    echo "samtools not found on PATH!"
    exit 1
fi

samtools idxstats ${bam} \
| awk -v bam=$bam 'BEGIN{OFS="\t"; FS="\t"} {
    if($1 == "*"){
        unmapped=$3+$4
    } else {
        mapped+=$3
        mateun+=$4
        refsize+=$2
    }
}END{print bam, mapped, mateun, unmapped, refsize}'  
done

exit 0
