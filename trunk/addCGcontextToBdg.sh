#!/usr/bin/env bash

set -e
set -o pipefail

bdg=$1
ref=$2
strand=$3

if [ "$bdg" = "-h" ] || [ "$bdg" = "--help" ] || [ "$bdg" = "" ] || [ "$ref" = "" ] || [ "$strand" = ""  ]; then
echo "
DESCRIPTION
    Add (annotate) CG context to bedgraph file produced by bam2methylation.py. Input format
    is tab separated:
        chrom, start, end, pct_met, cnt_met, tot, strand
    What matters is:
    - 'chrom' and 'start' matches the C to annotate. Typically 'end' is start+1 but
        it doesn't have to
    
REQUIREMENTS
    - bedtools on path
    - fastaFromBed.py on path (bedtools fastaFromBed will do the same but very slow for
        long files)
    - Fasta file indexed with index named <seqname>.fai in the same dir as input fasta.
      Use e.g 'samtools faidx' to create it.

USAGE
    addCGcontextToBdg.sh <bedgraph> <fasta-ref> <strand-column-index>

NOTES:
    Error handling is not very robust, if one of the pipes fails the script keeps going!

Version: 0.2
"
exit 1
fi


awk -v strandCol=$strand 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", ".", $strandCol}' $bdg \
| awk '{if ($6 != "+" && $6 != "-") {print "Strand information must be coded as + and -" > "/dev/stderr"; exit 1} else {print $0}}' \
| slopBed -l 0 -r 2 -s -i - -g ${ref}.fai \
| fastaFromBed.py -bed - -fi $ref -fo - -tab -s \
| awk 'BEGIN{IGNORECASE=1} {if ($2~/^CG/)
                                {print "CG"}
                            else if ($2~/^C[ACTU][G]/)
                                {print "CHH"}
                            else if ($2~/^C[ACTU][ACTU]/)
                                {print "CHH"}
                            else if ($2~/^C/)
                                {print "CNN"}
                            else
                                {print "ERROR: There is no C at this position!" > "/dev/stderr"; exit 1}}' \
| paste $bdg -
