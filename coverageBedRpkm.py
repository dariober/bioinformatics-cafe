#!/usr/local/bin/python


import sys
import subprocess
import shlex

if len(sys.argv) == 1:
    sys.exit("""
Execute coverageBed (in default mode) and adds to the output an additional
column of rpkm value. The RPKM column is inserted *before* the read count column
(in this way the output of coverageBedRpkm.py has the same format os coverageBed)

NOTES:
    * The read count is taken from the *4th* last column.
    * coverageBed is assumed to be on the path.
    * The output from coverageBed is loaded in memory
    
USAGE:

    coverageBedRpkm.py '<args to coveregeBed>' > <outfile>

TODO:
    
    - Add argument to get read-count from coverageBed from columns other then
    the 4th last

EXAMPLE:
    ## Note the single quotes around the first arg:
    coverageBedRpkm.py '-abam aln.bam -b ref.bed' > outcoverage.bed    
    """)

if len(sys.argv) != 2:
    sys.exit("""%s: Exactly 1 positional argument is required. E.g.
%s '-abam aln.bam -b ref.bed' > outcoverage.bed
             """ %(sys.argv[0], sys.argv[0]))

" ------------------------------[ Constants ]--------------------------------- "

READCOL= -4 ## Where to get the read-count in the output of coverageBed (-4 means 4th from right)

" ------------------------[ Functions ]--------------------------------------- "

def get_rpkm(readcount, mappedreads, regionlength):
    """
    Calculate rpkm as in http://www.clcbio.com/manual/genomics/Definition_RPKM.html
    """
    if regionlength == 0:
        return(0)
    rpkm= float(mappedreads) / ((readcount/1000000.0) * (regionlength/1000.0))
    return(rpkm)

## Count reads in bam file
args= shlex.split(sys.argv[1])
try:
    abam_ix= args.index('-abam')
    abam= args[abam_ix + 1]
except ValueError:
    abam= None
try:
    a_ix= args.index('-a')
    a= args[a_ix + 1]
except ValueError:
    a= None

if a is None and abam is None:
    sys.exit('%s: I cannot find the -abam or -a argumnet in the string %s' %(sys.argv[0], sys.argv[1]))
elif a is not None and abam is not None:
    sys.exit('%s: Both the -abam or -a argumnets found in string %s' %(sys.argv[0], sys.argv[1]))
elif abam is not None:
    p= subprocess.Popen('samtools view -c %s' %(abam), shell= True, stdout=subprocess.PIPE)
    readcount= int((p.communicate()[0]).strip())
elif a is not None:
    p= subprocess.Popen("""wc %s | awk {'print $1'}""" %(a), shell= True, stdout=subprocess.PIPE)
    readcount= int((p.communicate()[0]).strip())
else:
    sys.exit('%s: Unexpected conditions' %(sys.argv[0]))

## Execute coverageBed
cmd= 'coverageBed ' + sys.argv[1]
p= subprocess.Popen(cmd, shell= True, stdout=subprocess.PIPE)
coverage= p.communicate()[0]
coverage= coverage.split('\n')

for line in coverage:
    if line == '':
        continue
    line= line.rstrip('\n\r')
    line= line.split('\t')
    regionlength= int(line[2]) - int(line[1])
    mappedreads= int(line[READCOL])
    rpkm= get_rpkm(readcount, mappedreads, regionlength)
    line.insert(READCOL, rpkm)
    line= [str(x) for x in line]
    line= '\t'.join(line)
    print(line)
sys.exit()
