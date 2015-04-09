#!/bin/bash

docstring="Peak merger and overlap report \\n
USAGE \\n
mergePeaks.sh <a.bed> [b.bed ...] \\n
\\n
REQUIRES on PATH \\n
- bedtools merge Version: v2.23.0+ \\n
- tableCat.py from /Users/berald01/svn_git/bioinformatics-cafe \\n
\\n
OUTPUT\n
- Merged bed to stdout \\n
- Regions unique and shared by each file to stderr \\n
"

beds="$*"

if [[ ${beds} = "" || ${beds} = "-h" || ${beds} == "--help" ]]
then
    echo -e $docstring
    exit 1
fi

mydir=$(mktemp -dt "$0")

tableCat.py -i $beds \
| awk -v OFS="\t" '{print $1, $2, $3, $NF}' \
| sort -k1,1 -k2,2n \
| mergeBed -c 4,4 -o distinct,count_distinct -i - > ${mydir}/merged.bed 

cat ${mydir}/merged.bed 

cut -f4 ${mydir}/merged.bed \
| sort \
| uniq -c >&2

rm -r $mydir

exit
