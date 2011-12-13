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

filename_pileup= '/exports/work/vet_roslin_nextgen/dario/samtools/output/20100409_RNAseq_LPS_tophat/20100409_RNAseq_LPS_sscrofa9.56.pileup'

filename_slice= '/exports/work/vet_roslin_nextgen/dario/miscellanea/20100629/20100409_RNAseq_LPS_sscrofa9.56.pileup.slice'

contig= '4'

start= 94788700
end= 95009371

# -----------------------------------------------------------------------------

t0= time.time()

pileup= open(filename_pileup)
slice= open(filename_slice, 'w')

slice.write('rname\tpos\tcov\n')

base_counter= start
counter= 0

print('\nStart reading pileup file...\n\n')

for line in pileup:
    line= line.split('\t')
    pileup_pos = int(line[1])
    if (line[0] == contig) and (pileup_pos >= start) and (pileup_pos <= end):
        while pileup_pos > base_counter:
            """ Add a line with 0 count if a position is not covered. I.e. has no entry in the pileup """ 
            missing_pos= [contig, str(base_counter), '0']
            missing_pos= '\t'.join(missing_pos)
            missing_pos= missing_pos + '\n'
            slice.write(missing_pos)
            base_counter += 1
        outline= [line[0], line[1], line[7]]
        outline= '\t'.join(outline)
        slice.write(outline + '\n')
        base_counter += 1
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



