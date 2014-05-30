#!/usr/bin/env bash

bdg=$1
ref=$2

if [ "$bdg" = "-h" ] || [ "$bdg" = "--help" ] || [ "$bdg" = "" ] || [ "$ref" = "" ]; then
echo "
DESCRIPTION
    Add (annotate) CG context to bedgraph file produced by bam2methylation.py. Input format
    is tab separated:
        chrom, start, end, pct_met, cnt_met, tot, strand
    What matters is:
    - 'chrom' and 'start' matches the C to annotate. Typically 'end' is start+1 but
        it doesn't have to
    - Strand information is in **7th column** NOT 6th as in standard bed (sorry).
    
REQUIREMENTS
    - bedtools on path
    - fastaFromBed.py on path (bedtools fastaFromBed will do the same but very slow for
        long files)
    - Fasta file indexed with index named `<seqname>.fai` in the same dir as input fasta.
      Use e.g `samtools faidx` to create it.

USAGE
    addCGcontextToBdg.sh <bedgraph> <fasta-ref> 
"
exit 1
fi

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $7, $6}' $bdg \
| slopBed -l 0 -r 2 -s -i - -g ${ref}.fai \
| fastaFromBed.py -q -bed - -fi $ref -fo - -tab -s \
| awk 'BEGIN{IGNORECASE=1} {if ($2~/^CG/)
                                {print "CG"}
                            else if ($2~/^C[ACTU][G]/)
                                {print "CHH"}
                            else if ($2~/^C[ACTU][ACTU]/)
                                {print "CHH"}
                            else if ($2~/^C../)
                                {print "CNN"}
                            else
                                {print "NNN"}}' \
| paste $bdg -