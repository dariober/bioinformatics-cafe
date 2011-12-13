"""
  Intersection between refernce transcriptome and assembled transcriptome
"""

import interval
    
f_exon_intervals= 'C:/Tritume/exon_intervals.txt'
exon_intervals= (open(f_exon_intervals, 'r')).readlines()
exon_intervals= [x.split('\t') for x in exon_intervals]
interval_assem= [eval(x[2]) for x in exon_intervals]
interval_ref= [eval(x[3]) for x in exon_intervals]

intersect_exon= []
for i in xrange(len(interval_assem)):
    intersect_exon.append(interval_assem[i] & interval_ref[i])
intersect_sizes= []
for int in intersect_exon:
    intersect_sizes.append([(x[1]-x[0]) for x in int])
overlap= [sum(x) for x in intersect_sizes]

ref_sizes=[]
for int in interval_ref:
    ref_sizes.append([(x[1]-x[0]) for x in int])
tot_ref_sizes= [sum(x) for x in ref_sizes]

assem_sizes=[]
for int in interval_assem:
    assem_sizes.append([(x[1]-x[0]) for x in int])
tot_assem_sizes= [sum(x) for x in assem_sizes]

perc_coverage= round((sum(overlap)/sum(tot_ref_sizes))*100, 2)
print('Total size of reference transcripts: ' + str(sum(tot_ref_sizes)))
print('Total size of assembled transcripts: ' + str(sum(tot_assem_sizes)))
print('Total size of overlap: ' + str(sum(overlap)) + ' (' + str(perc_coverage) + '% of the reference)')

"""

exon_intervals[2166]
assem_sizes[2166]
tot_assem_sizes[2166]

ref_sizes[2166]
tot_ref_sizes[2166]

"""