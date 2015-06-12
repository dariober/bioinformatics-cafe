#!/usr/bin/env bash

set -e
set -o pipefail

VERSION='0.1.0'

bdg=$1
ref=$2
strand=$3
l=$4
r=$5

if [ "$bdg" = "-h" ] || [ "$bdg" = "--help" ] || [ "$bdg" = "" ] || [ "$ref" = "" ] || [ "$strand" = ""  ] || [ "$l" = ""  ] || [ "$r" = ""  ]; then
echo "
DESCRIPTION
    Add sequence context to intervals in a bed file

POSITIONAL ARGUMENTS
1. Bed file to annotate 
2. Reference fasta file. *.fai index must be present in the same dir.
3. Column index (1-based) in bed file with strand info. '+' for forward '-' for reverse
4. Add this many bases to from UPstream of the interval (i.e. to the left if +strand on forw)
5. Add this many bases to from DOWNstream of the C (i.e. to the right if -strand on forw)

REQUIREMENTS
    - bedtools on path
    - fastaFromBed.py on path (bedtools fastaFromBed will do the same but very slow for
        long files)
    - Fasta file indexed with index named <seqname>.fai in the same dir as input fasta.
      Use e.g 'samtools faidx' to create it.

USAGE
    addContextToBed.sh <bedgraph> <fasta-ref> <strand-column-index> <up> <down>

NOTES:
    Error handling is not very robust, if one of the pipes fails the script keeps going!

Version: $VERSION
"
exit 1
fi

awk -v strandCol=$strand 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", ".", $strandCol}' $bdg \
| awk '{if ($6 != "+" && $6 != "-") {print "Strand information must be coded as + and -" > "/dev/stderr"; exit 1} else {print $0}}' \
| slopBed -l $l -r $r -s -i - -g ${ref}.fai \
| fastaFromBed.py -bed - -fi $ref -fo - -tab -s 2> /dev/null \
| cut -f 2 \
| paste $bdg -


