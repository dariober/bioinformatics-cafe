#!/usr/bin/env bash

VERSION='0.1.0'

if [[ $1 == '-h' || $1 == '--help' || $1 == '' ]]
then
    echo ""
    echo "Interleave one or more pairs of fastq files."
    echo ""
    echo "USAGE"
    echo "catInterleaveFastq.sh <R1.L001.fq.gz,R1.L002.fq.gz,...> <R2.L001.fq.gz,R2.L002.fq.gz,...>"
    echo ""
    echo "First and second arguments are comma separated lists of first-in-pair and "
    echo "second-in-pair fastq file(s), respectively. No spaces or commas are allowed "
    echo "within filenames. Files can be either ALL decompressed or ALL gzipped "
    echo "compressed. If gzipped they must have extension '.gz'."
    echo "There is no check whether files are correctly paired and well formed."
    echo ""
    echo "Version $VERSION"
    exit 0
fi

set -ef -o pipefail

if [[ $1 == "" || $2 == "" || $3 != "" ]]
then
    echo "Please provide 2 positional arguments"
    exit 1
fi

r1=`echo $1 | sed 's/,/ /g'`
r2=`echo $2 | sed 's/,/ /g'`

# Some checks
# -----------
for x in $r1 $r2
do
    if [[ ! -f $x ]]; 
    then
        echo "Error: file $x not found."
        exit 1
    fi
done

nr1=`echo $r1 | wc -w`
nr2=`echo $r2 | wc -w`
if [[ $nr1 != $nr2 ]]
then
    echo "Error: number of first-in-pair and second-in-pair files are not the same."
    exit 1
fi

# Start processing
# ----------------

if [[ $1 == *.gz && $2 == *.gz ]]
then
    paste <(gzip -c -d $r1 | paste - - - -) \
          <(gzip -c -d $r2 | paste - - - -) \
    | tr '\t' '\n'
else
    paste <(cat $r1 | paste - - - -) \
          <(cat $r2 | paste - - - -) \
    | tr '\t' '\n'
fi
