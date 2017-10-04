#!/usr/bin/env bash

VERSION='0.1.0'

if [[ $1 == '-h' || $1 == '--help' || $1 == '' ]]
then
    echo ""
    echo "Interleave one or more pairs of fastq files."
    echo ""
    echo "USAGE"
    echo "catInterleaveFastq.sh -1 R1.L001.fq.gz R1.L002.fq.gz ... -2 R2.L001.fq.gz R2.L002.fq.gz ..."
    echo ""
    echo "First and second arguments are comma separated lists of first-in-pair and "
    echo "second-in-pair fastq file(s), respectively. No spaces are allowed "
    echo "within filenames. Files can be either ALL decompressed or ALL gzipped "
    echo "compressed. If gzipped they must have extension '.gz'."
    echo "There is no check whether files are correctly paired and well formed."
    echo ""
    echo "Version $VERSION"
    exit 0
fi

set -euf -o pipefail

# Collect args
r1=""
r2=""
xstart=0
for x in $@
do
    if [[ $x == '-1' ]]
    then
        xstart=1 # Start collecting 1st-in-pair files
        continue
    fi
    
    if [[ $x == '-2' ]]
    then
        xstart=2 # Start collecting 2nd-in-pair files
        continue
    fi

    if [[ ! -f $x ]]; 
    then
        echo "Error: file $x not found."
        exit 1
    fi

    if [[ $xstart == 1 ]]
    then
        r1="$r1 $x"
    elif [[ $xstart == 2 ]]
    then
        r2="$r2 $x"
    fi
done

if [[ $r1 == "" || $r2 == "" ]]
then
    echo "Please provide arguments to parameters -1 and -2"
    exit 1
fi

# Some checks
# -----------
nr1=`echo $r1 | wc -w`
nr2=`echo $r2 | wc -w`
if [[ $nr1 != $nr2 ]]
then
    echo "Error: number of first-in-pair and second-in-pair files are not the same."
    exit 1
fi

# Start processing
# ----------------

if [[ $r1 == *.gz && $r2 == *.gz ]]
then
    paste <(gzip -c -d $r1 | paste - - - -) \
          <(gzip -c -d $r2 | paste - - - -) \
    | tr '\t' '\n'
else
    paste <(cat $r1 | paste - - - -) \
          <(cat $r2 | paste - - - -) \
    | tr '\t' '\n'
fi
