#! /usr/bin/python

"""
  
  From a file in pileup format (samtools), extract the lines (genome positions) between 
  to coordinates (on the same contig/chromosome).

  The pileup should like like:

  MT      1       C       C       30      0       60      1       ^~.     a       ~
  MT      2       A       A       48      0       60      7       .^~,^~,^~,^~,^~,^~,     aX^]UV\ ~~~~~~~

  With:
  col[0]= rname (contig name)
  col[1]= position
  col[7]= coverage (1 and 7 in the sample above)

"""

import time 

# ------------------------------[ User input ]---------------------------------

filename_pileup= '/exports/work/vet_roslin_nextgen/dario/samtools/output/20100630_RNAseq_CTRL/20100630_RNAseq_CTRL.pileup'

filename_slice= '/exports/work/vet_roslin_nextgen/dario/miscellanea/20100701/20100630_RNAseq_CTRL_chr4.pileup'

contig= '4'

# -----------------------------------------------------------------------------

t0= time.time()

pileup= open(filename_pileup)
slice= open(filename_slice, 'w')

slice.write('rname\tpos\tref_base\t_read_base\tphred_consensus\tphred_snp\trms_quality\tno_reads\n')

counter= 0
counter_out= 0
print('\nStart reading pileup file...\n\n')

for line in pileup:
    line= line.split('\t')
    if line[0] == contig:
        outline= line[0:8]
        outline= '\t'.join(outline)
        slice.write(outline+'\n')
        counter_out += 1
    counter += 1
    if counter % 1000000 == 0:
        print('Read pilup file number: ' + str(counter))
#    if counter > 2000:
#        break
pileup.close()
slice.close()

t1= time.time()

print('\n' + 'Finished extracting slice in ' + str(round(t1-t0, 2)) + ' sec' + '\n')
print('Number of lines read in pileup input: ' + str(counter) + '\n')
print('Number of lines read in output file: ' + str(counter_out) + '\n')