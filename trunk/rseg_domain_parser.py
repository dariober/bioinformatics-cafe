#!/usr/local/bin/python

import sys
import re
import os

if len(sys.argv) != 2 or sys.argv[1] in ['-h', '--help', '-help']:
    sys.exit("""
%(fname)s:

DESCRIPTION 
    Extracts the ENRICHED regions from the output file "*.domains.bed"
    prooduced by rseg/rseg-diff. The extracted lines go to file *.enriched.bed.

    Prints to stdout some statistics relative to the enriched regions.

USAGE
    %(fname)s <rsegout-domains.bed>

    """ %{'fname': os.path.split(sys.argv[0])[1]})


def median(s):
    """
    Calculate the median of a list of numbers
    Found at http://bytes.com/topic/python/answers/852576-how-do-you-find-median-list
    """
    s.sort()
    i = len(s)
    if not i%2:
        return( (s[(i/2)-1]+s[i/2])/2.0)
    return(s[i/2])

fin= sys.argv[1]
fout= re.sub('domains\.bed$', 'enriched.bed', fin)

fhin= open(fin)
fhout= open(fout, 'w')

# ----------------------- [ Stuff to store ]-----------------------------------

domain_states= {'ENRICHED': 0}         ## Dictionary to hold the count of states {'ENRICHED': x, 'BACKGROND': y, 'UNCERTAIN': z, ...}
enriched_sizes= []        ## List for the size of the enriched regions
enriched_counts= []    ## 5th column from domain file. From manuka: The 5th column gives the average read count in the domain
enriched_domain_score= [] ## 6th column: The sum of posterior scores of all bins within this domain; it measures both the quality and size of the domain

for line in fhin:
    line= line.rstrip('\n\r').split('\t')
    domain_states[line[3]]= domain_states.get(line[3], 0) + 1
    if line[3] == 'ENRICHED':
        enriched_sizes.append(int(line[2]) - int(line[1]))
        enriched_counts.append(float(line[4]))
        enriched_domain_score.append(float(line[5]))
        fhout.write('\t'.join(line) + '\n')
fhin.close()
fhout.close()

# -------------------------[ Stats ]-------------------------------------------

domain_counts= []
for k in sorted(domain_states.keys()):
    c= k + ':' + str(domain_states[k])
    domain_counts.append(c)
domain_counts= '; '.join(domain_counts)

if domain_states['ENRICHED'] > 0:
    avg_enriched_size= str(round(sum(enriched_sizes) / float(len(enriched_sizes)), 0))
    med_enriched_size= str(median(enriched_sizes))
    avg_enriched_counts= str(round(sum(enriched_counts) / float(len(enriched_counts)), 0))
    avg_enriched_domain_score= str(round(sum(enriched_domain_score) / float(len(enriched_domain_score)), 0))
    enriched_length= str(sum(enriched_sizes))
else:
    avg_enriched_size= 'NA'
    med_enriched_size= 'NA'
    avg_enriched_counts= 'NA'
    avg_enriched_domain_score= 'NA'
    enriched_length= 'NA'

outline= [sys.argv[1], domain_counts, 'AVG_ENRICHED_SIZE', avg_enriched_size, 'MEDIAN_ENRICHED_SIZE', med_enriched_size, 'ENRICHED_LENGTH', enriched_length ,'AVG_COUNTS', avg_enriched_counts, 'AVG_SCORE', avg_enriched_domain_score]
print('\t'.join(outline))
sys.exit()
