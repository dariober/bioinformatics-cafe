#!/usr/bin/env python

import sys
import subprocess
import time
import argparse

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Submits jobs to lsf in batches so that no more than n jobs run at the same time

USAGE
    LSFbatchSubmit.py <N> [<job.1> <job.2> ... <job.x> | -]

N: Number of concurrent jobs
job.1 ... : Script files that will be submitted. Use - to read the list of files from stdin

EXAMPLE:
    LSFbatchSubmit.py -j test- -n 4 -f H1.GSM882245.SRR449051.fastq.gz.sh ...


""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--files', '-f',
                   required= True,
                   nargs= '+',
                   help='''List of scripts to submit to LSF. Use - to read the list
from stdin.
                   ''')

parser.add_argument('--njobs', '-n',
                   required= False,
                   default= 1,
                   type= int,
                   help='''Run up to this many jobs at the same time. Default 1'''
                   )

parser.add_argument('--jobpre', '-j',
                   required= True,
                   help='''A prefix to prepend to each file name. The resulting str
will be used as job identifier and to name the bsub log file. E.g.:
LSFbatchSubmit.py --files script.sh --jobpre bsj- >>> bsub -J bsj-script.sh -oo bsjb-script.sh.bsub.log ...

MEMO: LSFbatchSubmit.py will use jobpre to get the number of running jobs using
bjobs -J "jobpre*"'''
                   )

parser.add_argument('--bsubOpt',
                   required= False,
                   default= '',
                   help='''String of further arguments to pass to bsub. Do not include
-J and -oo as these are automatically generated from --jobpre. Consider `-R "rusage[mem=xxxx]"`.
Default is '' (no further options)'''
                   )

args= parser.parse_args()

# -----------------------------------------------------------------------------
def getJobs(x):
    """ Get the number of jobs matching the pattern x"""
    cmd= 'bjobs -w -J %s' %(x)
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    p.wait()
    j= p.stdout.readlines()
    e= p.stderr.read().strip()
    if not e == '':
        sys.exit('I cannot get a list of running jobs')
    return(j)
# -----------------------------------------------------------------------------

if args.njobs <= 0:
    sys.exit('Number of jobs must be an integer > 0')

if args.files == ['-']:
    jobs= sys.stdin.readlines()
    jobs= [x.strip() for x in jobs]
else:
    jobs= args.files

for f in jobs:
    try:
        x= open(f)
    except:
        sys.exit('File %s cannot be opened.' %(f))
    x.close()

## Start submitting jobs
for f in jobs:
    cmd= "bsub -J %(jobid)s%(sh)s -oo %(sh)s.bsub.log %(bsubOpt)s < %(sh)s" %{'jobid': args.jobpre, 'sh': f, 'bsubOpt': args.bsubOpt}
    print(cmd)
    p= subprocess.Popen(cmd, shell= True)
    p.wait()
#    for j in subjobs:
#        print(j.strip())
    subjobs= getJobs(args.jobpre + '*')
    while len(subjobs) > args.njobs: ## MEMO: 1st lie of subjobs is the header
        time.sleep(60)
        subjobs= getJobs(args.jobpre + '*')
sys.exit()



