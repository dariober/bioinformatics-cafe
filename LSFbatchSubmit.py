#!/usr/bin/env python

import sys
import subprocess
import time
import argparse
import glob
import re

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Submits jobs to lsf in batches so that no more than n jobs run at the same time.
    Say you have several sh scripts to submit to lsf which could be run in parallel.
    If however you can't run them in parallel because, for example, there is not enough disk space,
    use LSFbatchSubmit to run at most N jobs (scripts) at the same time.
    
    Jobs are submitted in the same order as they are supplied. Therefore, with n=1
    a pipeline of jobs can be executed (sort of the same as using 'bsub -w "ended(...)"').

EXAMPLE:
    ## Execute all scripts matching H1.4*.fastq.gz.sh no more than two at a time
    LSFbatchSubmit.py -j bismark_pipeline- -n 2 -f H1.4*.fastq.gz.sh --echo
        bsub -J bismark_pipeline-H1.45.fastq.gz.sh -oo H1.45.fastq.gz.sh.bsub.log  < H1.45.fastq.gz.sh
        bsub -J bismark_pipeline-H1.46.fastq.gz.sh -oo H1.46.fastq.gz.sh.bsub.log  < H1.46.fastq.gz.sh
        bsub -J bismark_pipeline-H1.47.fastq.gz.sh -oo H1.47.fastq.gz.sh.bsub.log  < H1.47.fastq.gz.sh
        bsub -J bismark_pipeline-H1.48.fastq.gz.sh -oo H1.48.fastq.gz.sh.bsub.log  < H1.48.fastq.gz.sh


""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--files', '-f',
                   required= True,
                   nargs= '+',
                   help='''List of scripts to submit to LSF. Use - to read the list
from stdin. Wild cards in file names are allowed and are expanded using glob.glob().
Duplicate files e.g. from '-f script.1.sh script.*.sh' are allowed and only one file
will be submitted.
If the file list contains non-existant files, LSFbatch will exit with error.
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

parser.add_argument('--echo', '-e',
                    action= 'store_true',
                    help= '''Do not submit jobs. Only print the command that would be executed ('bsub -J ... -oo ...')
                    ''')

parser.add_argument('--time', '-t',
                    required= False,
                    type= int,
                    default= 60,
                    help= '''Check the job list by means of bjobs every this many seconds (default 60).
                    ''')

args= parser.parse_args()

# -----------------------------------------------------------------------------
def getJobs(x):
    """ Get the number of jobs matching the pattern x
    Returns a list of strings, one for each line of bjobs output or empty list []
    if the pattern doesn't match anything.
    
    Output for matching patterns looks like:
    ['JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME',
     '85894   berald01 RUN   cluster    uk-cri-lcst01 crinode52   bismark_pipeline-H1.GSM882245.SRR449052.fastq.gz.sh Jan  4 01:03',
     '87622   berald01 RUN   cluster    uk-cri-lcst01 crinode52   bismark_pipeline-H1.GSM882245.SRR449055.fastq.gz.sh Jan  4 09:02',
     ...]

    """
    cmd= 'bjobs -w -J %s' %(x)
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    p.wait()
    j= p.stdout.readlines()
    j= [x.strip() for x in j]
    e= p.stderr.read().strip()
    if not e == '':
        if re.match('Job <.*> is not found', e) or e == 'No unfinished jobs found.':
            return([])
        else:
            print('Error executing %s' %(cmd))
            print(e)
            sys.exit(1)
    return(j)
# -----------------------------------------------------------------------------

if args.njobs <= 0:
    print('Number of jobs must be an integer > 0')
    sys.exit(1)

if args.files == ['-']:
    jobs= sys.stdin.readlines()
    jobs= [x.strip() for x in jobs]
else:
    jobs= []
    for f in args.files:
        globf= glob.glob(f)
        if globf == []:
            sys.exit('Cannot find file %s' %(f))
        jobs= jobs + globf
## Remove duplicates while preserving the order.
jobsUnique= []
for j in jobs:
    if j not in jobsUnique:
        jobsUnique.append(j)
jobs= [x for x in jobsUnique]

## Start submitting jobs
for f in jobs:
    cmd= "bsub -J %(jobid)s%(sh)s -oo %(sh)s.bsub.log %(bsubOpt)s < %(sh)s" %{'jobid': args.jobpre, 'sh': f, 'bsubOpt': args.bsubOpt}
    print(cmd)
    if args.echo:
        continue
    subjobs= getJobs(args.jobpre + '*')
    while len(subjobs) > args.njobs: ## MEMO: 1st line of subjobs is the header
        time.sleep(args.time)
        subjobs= getJobs(args.jobpre + '*')
    p= subprocess.Popen(cmd, shell= True)
    p.wait()
    if f == jobs[-1]:
        ## Exit loop if the last job has been submitted
        break
sys.exit()



