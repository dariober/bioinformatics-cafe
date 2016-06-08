#!/bin/bash

VERSION='0.1.0'

docstring="\n
Extract genome file from bam. Useful for tools that require a genome file \n
as input, like intersectBed \n
\n
USAGE \n
genomeFileFromBam.sh <in.bam> \n
\n
REQUIRES: samtools and Unix tools: awk, sed, paste \n

version: $VERSION \n
"

inbam=$1

if [ "$inbam" = "-h" ] || [ "$inbam" = "--help" ] || [ "$inbam" = "" ]; then
    echo -e $docstring
    exit 1
fi

paste \
<(samtools view -H $inbam \
| awk '$0 ~ "^@SQ"' \
| sed -r 's/^@SQ.+?SN://' \
| sed 's/\t.*//') \
<(samtools view -H $inbam \
| awk '$0 ~ "^@SQ"' \
| sed -r 's/^@SQ.+?LN://' \
| sed 's/\t.*//')