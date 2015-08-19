#!/bin/bash

set -e
set -o pipefail

docstring='
\n
Parse output of samtools idxstats to get alignment stats from one or more bam files.
\n
Requires awk and samtools on PATH.
\n\n
USAGE \n
readCountFromBam.sh <aln.bam> <aln2.bam> ... \n

\n
OUTPUT \n
Tab separated table with columns:\n
1. Input file name \n
2. Count mapped. Sum 3rd col [samtools view -F4 -c aln.bam] \n
3. Count reads unmapped, mate mapped. Sum 4th col [samtools view -f4 -F8 -c aln.bam]  \n
4. Count unmapped (\* line) [samtools view -f4 -f8 -c aln.bam] \n
5. Ref size (sum 2nd col) \n
6. Percent mapped. 100 x col#2 / (col#2 + col#3 + col#4)
\n\n
Note that a read unmapped but with mate mapped (col #3) is assigned to the \n
same position as the mate (at last by bwa). This means that a read can be \n
flagged as unmapped while having chrom and pos assigned. \n\n

Version 0.2.0'

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
}END{print bam, mapped, mateun, unmapped, refsize, 100 * mapped / (mapped + mateun + unmapped)}'  
done

exit 0
