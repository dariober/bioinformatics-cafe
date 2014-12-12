#!/bin/bash

set -e
set -o pipefail

docstring="Get read count from BAM file.
\n\n
This script parses the output of samtools idxstats to quickly get the total read
count from a bam file.
\n\n
Requires awk and samtools on PATH.
\n\n
USAGE \n
readCountFromBam.sh <aln.bam>
\n\n
OUTPUT \n
* Input file name
# mapped (sum 3rd col) \n
# mate unmapped (sum 4th col) \n
# unmapped (* line) \n
# ref size (sum 2nd col) \n
"

bam=$1

if [[ ${bam} = "" || ${bam} = "-h" || ${bam} == "--help" ]]
then
    echo -e $docstring
    exit 1
fi

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

exit 0
