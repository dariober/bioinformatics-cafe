#!/cygdrive/c/Python26/python

"""
Given a list of pileup files, merge the files to have
a sort of crosstab format with row headers rname, pos, ref_base
and column headers the read counts from each file. Positions
present absent in one pileup are reported as 0:

rname  pos  ref_base  pile1  pile2 ... pileN
1      1    A         10     4         0
1      2    C         0      2         4
...

Edit function parse_pileup() to have different columns in output.

"""
from flatten import flatten
from merge_join import merge_join
import os
import sys

# --------------------------[ Input/Output ]-----------------------------------

# Set working dir:
os.chdir('F:/data/20100630_RNAseq_pileup/')

# List of files to pileup. These files, parsed, will be held in memory.
pileup_f= ['20100630_RNAseq_CTRL_contig.pileup/20100630_RNAseq_CTRL_1.pileup',
           '20100630_RNAseq_LPS_contig.pileup/20100630_RNAseq_LPS_1.pileup']

# Does the pileup files have header line(True/False)? If True, it will be discarded:
header= True
           
# List of names to use as column headers (row headers plus one for each file above):
colnames= ['rname', 'pos', 'ref_base', 'n_ctrl', 'n_lps']

# Output file:
out_pileup= '20100630_RNAseq_1.merge.pileup'

# -----------------------------------------------------------------------------

def parse_pileup(pileup):
    """ Takes an a pileup file (E.g. pileup= 'myfile.pileup') and returns a list 
    of key/values suitable for merge_join. They have to look like:
    [(('X', 1, 'A'), 15), (...)] meaning [((chromosome, position, ref_base), count), (...)]
    """
    pilefile= open(pileup)
    if header == True:
        pilefile.readline() # Discard header line
    t_pileup= []        # Where the list of key/values will be stored
    for line in pilefile:
        line= line.rstrip('\n').split('\t')
        key= (line[0], int(line[1]), line[2])
        value= int(line[7])
        t_pileup.append((key, value))
    pilefile.close()
    return(t_pileup) 
ct_pileup= open(out_pileup, 'w')
ct_pileup.write('\t'.join(colnames).rstrip('\n') + '\n')

# Change this manually to add more :
print('Parsing file ' + pileup_f[0])
p1= parse_pileup(pileup_f[0])
print('Parsing file ' + pileup_f[1])
p2= parse_pileup(pileup_f[1])

# And this add as many p1, p2, ... pN as necessary:
n= 0
for t in merge_join(p1, p2):
    t= list(flatten(t))
    for i in range(0, len(t)):
        if t[i] == None:
            t[i]= 0
    t= [str(i) for i in t]
    ct_pileup.write('\t'.join(t) + '\n')
    if n % 250000 == 0:
        print(str(n) + ' lines merged.')
ct_pileup.close()

 

