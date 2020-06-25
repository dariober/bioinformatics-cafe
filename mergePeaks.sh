#!/bin/bash

docstring="Peak merger and overlap report \\n
USAGE & EXAMPLE\\n
mergePeaks.sh <a.bed> [b.bed ...] \\n
\\n
mergePeaks.sh a.bed b.bed c.bed > merged.bed 2> summary.txt
\\n
\\n
REQUIRES on PATH \\n
- bedtools merge Version: v2.23.0+ \\n
- tableCat.py from https://github.com/dariober/bioinformatics-cafe \\n
\\n
OUTPUT\n
- Merged bed to stdout \\n
- Regions unique and shared by each file to stderr \\n

Version 0.3.0
"

beds="$*"

if [[ ${beds} = "" || ${beds} = "-h" || ${beds} == "--help" ]]
then
    echo -e $docstring
    exit 1
fi

tmp=$(mktemp tmp.XXXXXXXXXX)

for bed in $beds
do
    awk -v OFS='\t' '{print $0, FILENAME}' $bed
done \
| grep -v '^#' \
| awk -v OFS="\t" '{print $1, $2, $3, $NF}' \
| sort -k1,1 -k2,2n \
| mergeBed -c 4,4 -o distinct,count_distinct -i - > ${tmp}

cat ${tmp}

cut -f4 ${tmp} \
| sort \
| uniq -c >&2

rm $tmp

exit
