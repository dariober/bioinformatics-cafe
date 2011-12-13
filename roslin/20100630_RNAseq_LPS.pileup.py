#! /usr/bin/python

"""
  
  From a file in pileup format (samtools), extracts the individual contigs (chromosomes).

  The pileup should like like:

  MT      1       C       C       30      0       60      1       ^~.     a       ~
  MT      2       A       A       48      0       60      7       .^~,^~,^~,^~,^~,^~,     aX^]UV\ ~~~~~~~

  With:
  col[0]= rname (contig name)
  col[1]= position
  col[7]= coverage (1 and 7 in the sample above)

  The first 8 columns are extracted.

"""

import time 

# ------------------------------[ User input ]---------------------------------

filename_pileup= '/exports/work/vet_roslin_nextgen/dario/samtools/output/20100630_RNAseq_CTRL/20100630_RNAseq_CTRL.pileup'



## Path and base name for the extracted files. '.pileup' will be added as extension
filename_slice= '/exports/work/vet_roslin_nextgen/dario/miscellanea/20100705/20100630_RNAseq_CTRL_'

# -----------------------------------------------------------------------------


t0= time.time()

pileup= open(filename_pileup)


counter= 0
print('\nStart reading pileup file...\n\n')

line= pileup.readline()
line= line.split('\t')
while True:
    if line == '':
        break
    current_contig= line[0]
    current_outfile= open(filename_slice + current_contig + '.pileup', 'w')    
    current_outfile.write('rname\tpos\tref_base\t_read_base\tphred_consensus\tphred_snp\trms_quality\tno_reads\n')
    print('Extracting contig: ' + current_contig)
    while line[0] == current_contig:
        outline= line[0:8]
        outline= '\t'.join(outline)
        current_outfile.write(outline+'\n')
        line= pileup.readline()
        if line == '':
            break
        line= line.split('\t')
        counter += 1
    print('Positions in ' + current_contig + ': ' + str(counter) + '\n')
    conter= 0

pileup.close()
current_outfile.close()

t1= time.time()

print('\n' + 'Finished extracting contigs in ' + str(round(t1-t0, 2)) + ' sec' + '\n')
