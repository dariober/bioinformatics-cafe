#!/usr/bin/env python

import sys

VERSION='0.1.0'

docstring= """DESCRIPTION
Convert FIMO txt output to bed format with global coordinates. It assumes the
sequence name column in fimo is in the format chrom:start-end as produced for
example by 'bedtools getfasta'

USAGE
fimoToGlobalBed.py <fimo.txt | - for stdin>

Typical use case: You have chipseq regions to pass to fimo to see where binding 
motifs are. For this you use fastaFromBed to give to fimo. Fimo gives you the
coordinates relative to the inpuyt sequence, but you want to know where the 
binding motifs are in genomic coordinates.

*** Example Input:
#pattern name   sequence name   start   stop    strand  score   p-value q-value matched sequence
MA0139.1        chr9:20635125-20636125  500     518     -       27.5246 1.37e-12        7.27e-05        TGGCCACCAGGGGGCGCTA
MA0139.1        chr3:19593586-19594586  474     492     -       27.2131 3.61e-12        7.27e-05        TGGCCACCAGGGGGCGCTG
MA0139.1        chr5:28420189-28421189  480     498     -       27.2131 3.61e-12        7.27e-05        TGGCCACCAGGGGGCGCTG
MA0139.1        chr7:127813362-127814362        486     504     +       27.2131 3.61e-12        7.27e-05        TGGCCACCAGGGGGCGCTG

*** Output (note header removed)
chr9	20635624	20635643	MA0139.1	chr9:20635125-20636125	500	518	-	27.5246	1.37e-12	7.27e-05	TGGCCACCAGGGGGCGCTA
chr3	19594059	19594078	MA0139.1	chr3:19593586-19594586	474	492	-	27.2131	3.61e-12	7.27e-05	TGGCCACCAGGGGGCGCTG
chr5	28420668	28420687	MA0139.1	chr5:28420189-28421189	480	498	-	27.2131	3.61e-12	7.27e-05	TGGCCACCAGGGGGCGCTG
chr7	127813847	127813866	MA0139.1	chr7:127813362-127814362	486	504	+	27.2131	3.61e-12	7.27e-05	TGGCCACCAGGGGGCGCTG

Version %s
""" %(VERSION)

def seqnameToCoords(seqname):
    x= seqname[::-1].replace(':', '\t', 1).replace('-', '\t', 1)[::-1] # Replace last occurence of ':'  and '-' (chr name might contain these chars)
    lst= x.split('\t')
    assert len(lst) == 3
    lst[1]= int(lst[1])
    lst[2]= int(lst[2])
    return lst

if len(sys.argv) != 2 or sys.argv[1] in ('-h', '--help'):
    print docstring
    sys.exit(1)

if sys.argv[1] == '-':
    fin= sys.stdin
else:
    fin= open(sys.argv[1])


for line in fin:
    if line.strip().startswith('#'):
        continue
    line= line.strip().split('\t')
    locCoords= seqnameToCoords(line[1])
    globCoords= [locCoords[0], (locCoords[1]-1) + int(line[2]), locCoords[1] + int(line[3])]

    # A (weak) check that coords are correct: Span is the same size of the sequence
    assert (globCoords[2]-globCoords[1]) == len(line[8])

    globCoords= [str(x) for x in globCoords]
    outline= '\t'.join(globCoords + line)
    print outline
    
fin.close()
