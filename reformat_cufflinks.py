#!/home/berald01/.local/bin/python

"""
Reformat the output of cufflinks/transcripts.gtf to:
- Extract only 'transcript' rows
- Replace 'Cufflinks' in second column with an identifier
- Split attribute string to columns
- Add 'project_id' field
"""

import os
import sys
import re

id= sys.argv[1] ## 'ds0${i}'
infile= os.path.join('ds0' + id, 'transcripts.gtf')
project_id= '20120622_rnaseq_pdsa'
prefix= 'ds0'

if int(id) in [19,20,21]:     ## Make sure you assign the correct suffix!!
    suffix= '_ctrl_24h'
elif int(id) in [22,23,24]:
    suffix= '_pdsa_24h'
else:
    sys.exit('Unexpected index')
fin= open(infile)
for line in fin:
    line= line.strip().split('\t')
    if line[2] != 'transcript':
        continue
    else:
        line[1]= prefix + id + suffix ## Add sample_id
        idx= 8
        attr_col= line[8]
        for x in ['gene_id', 'transcript_id', 'FPKM', 'frac', 'conf_lo', 'conf_hi', 'cov']:
             p= x + ' "' + '(.*?)";'
             m= re.findall(p, attr_col)[0]
             line.insert(idx, m)
             idx += 1 
        line.append(project_id) ## Add project_id
        print('\t'.join(line))
fin.close()
sys.exit()