#!/usr/bin/env python

import sys
import pybedtools
import tempfile
import os

docstring= """Get the gene lengths from GTF file.
USAGE:
    geneLengthFromGTF.py genes.gtf

The GTF is expected to have attribute "gene_name". The gene length is the sum of the length
of the merged exons within "gene_name". Gene are returned in alphanumeric order.
"""

if len(sys.argv) != 2 or sys.argv[1] in ('-h', '--help'):
    sys.exit(docstring)

exons= pybedtools.BedTool(sys.argv[1]).filter(lambda x: x[2] == 'exon')
tmpexons= tempfile.NamedTemporaryFile(delete= False)
for line in exons:
    tmpexons.write('\t'.join([line.attrs['gene_name'], str(line.start), str(line.end)]) + '\n')
tmpexons.close()

merged= pybedtools.BedTool(tmpexons.name).sort().merge()

gene_length= {}
for line in merged:
    gene_length[line.chrom]= gene_length.get(line.chrom, 0) + line.length

for k in sorted(gene_length.keys()):
    print('\t'.join([k, str(gene_length[k])]))
sys.exit()