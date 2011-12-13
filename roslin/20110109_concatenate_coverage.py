#!/usr/bin/python

"""
Concatenate coverage files produced by BEDtools coverageBed. 
Concatenate only lines with non zero coverage (to reduce size)
Prepend to each line an identifier of the original file/library.
"""

## files to concatenate:
covfiles= ['LN_109.cov', 'LN_114.cov', 'LN_124.cov', 'LN_130.cov', 'LN_146.cov', 'LN_173.cov', 'LN_183.cov', 'LN_20B.cov', 'LN_21.cov', 'LN_38.cov', 'LN_47.cov', 'LN_50.cov', 'LN_57.cov', 'LN_58.cov', 'LN_92.cov']

## Output file:
covcat= open('LN_all.nonzero.cov', 'w')

## Prepend to each line the file name without .cov

# -----------------------------------------------------------------------------

for cov in covfiles:
    print('Start processing: ' + cov)
    lib_id= cov.rstrip('.cov')
    fcov= open(cov)
    for line in fcov:
        line= line.split('\t')
        if line[9] == '0':
            continue
        else:
            line= [lib_id] + line
            line= '\t'.join(line)
            covcat.write(line)
    fcov.close()
covcat.close()



