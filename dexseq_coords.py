#!/home/berald01/.local/bin/python

docstr= """
Read the gff output from dexseq_count.py to produce a table of exon IDs (gene_exon) and their
coordinates.

Example of gff attribute expected:
transcripts "NR_046018_1"; exonic_part_number "001"; gene_id "chr1_DDX11L1+"

Returned gene_exon ID: chr1_DDX11L1+:001

USAGE
    dexseq_coords.py <genes.gff>

SAMPLE OUTPUT:
chr1    11874   12227   chr1_DDX11L1+:001       +
chr1    12613   12721   chr1_DDX11L1+:002       +
chr1    13221   14408   chr1_DDX11L1+:003       +
chr1    14362   14829   chr1_WASH7P-:001        -
chr1    14970   15038   chr1_WASH7P-:002        -

"""

import sys
import re

def make_gene_exon_id(gff_attr):
    """Produce a gene_exon from the gff attribute
    """
    gff_attr= gff_attr.split('; ')
    gene_id= None
    exon_no= None
    for x in gff_attr:
        x= x.strip('"')
        if x.startswith('exonic_part_number'):
            exon_no= re.sub('exonic_part_number "', '', x)
        if x.startswith('gene_id'):
            gene_id= re.sub('gene_id "', '', x)
    if gene_id is None or exon_no is None:
        sys.exit('I cannot compose a valid gene_exon ID from %s.' %(gff_attr))
    gene_exon= gene_id + ':' + exon_no
    return(gene_exon)
    
if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
    print(docstr)
    sys.exit()

print('\t'.join(['chrom', 'start', 'end', 'gene_exon', 'strand']))
gff= open(sys.argv[1])
for line in gff:
    line= line.strip().split('\t')
    coord= []
    if line[2] == 'exonic_part':
        coord.append(line[0])
        coord.append(line[3])
        coord.append(line[4])
        coord.append(make_gene_exon_id(line[-1]))
        coord.append(line[6])
        print('\t'.join(coord))
gff.close()
sys.exit()