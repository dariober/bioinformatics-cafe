#!/usr/bin/env python3

import sys
import subprocess
import csv
import collections
from io import StringIO

docstring= """
sacct_easy.py [-tsv] [sacct-args]

Execute slurm sacct and return a friendly tabular output. 
Values in MaxRSS and ReqMem are in MB units.

OPTIONS

-tsv       
    Print output as tab separated values (better for further processing)

sacct-args 
    Further args to sacct, e.g. `-S 2017-09-10`. Do not include `--format`
"""

if '-h' in sys.argv or '--help' in sys.argv:
    print(docstring)
    sys.exit()

asTsv= False
if '-tsv' in sys.argv:
    asTsv= True
    sys.argv.remove('-tsv')

def normalizeMem(x):
    if x.strip() == '':
        return '';  
    x= x.strip('n')
    return round(int(x.replace("K", "000").replace("M", "000000").replace("G", "000000000"))/1e6)

def maxNameLen(sacct):
    xread= csv.DictReader(sacct.decode().split('\n'), delimiter='|')
    mlen= 10
    for line in xread:
        if len(line['JobName']) > mlen:
            mlen= len(line['JobName'])
    return mlen

def fillInJob(jobid, batchid):
    jobid['MaxRSS']= batchid['MaxRSS']
    return jobid    

def tabulate(line, maxlen, asTsv):
    if asTsv:
        return '\t'.join([str(x) for x in line])
    row_fmt= "{:<8}{:<%s}{:<12}{:<8}{:<8}{:<10}{:<10}{:<30}" % (maxlen + 2)
    return row_fmt.format(*line)

sacct= subprocess.check_output(['sacct', '--parsable2', '--format=JobID,JobName%50,Partition,MaxRSS,ReqMem,AllocCPUS,Elapsed,State'] + sys.argv[1:])

reader = csv.DictReader(sacct.decode().split('\n'), delimiter='|')

maxlen= maxNameLen(sacct)
print(tabulate(reader.fieldnames, maxlen, asTsv))

idjob= None
for line in reader:
    line['MaxRSS']= normalizeMem(line['MaxRSS'])
    line['ReqMem']= normalizeMem(line['ReqMem'])
    if '.batch' not in line['JobID']:
            if idjob is not None:
                # This line is not a batch job, so the previous line
                # can be printed as it doesn't have an associated batch job.
                lst= list(idjob.values())
                print(tabulate(lst, maxlen, asTsv))
            idjob= line
    else:
        # This is a batch job. So associate to it the job ID line    
        if idjob['JobID'] == line['JobID'].replace('.batch', ''):
            lst= list(fillInJob(idjob, line).values())
            print(tabulate(lst, maxlen, asTsv))
        else:
            print(idjob)
            print(line)
            raise Exception("Cannot find job id")
        idjob= None

if idjob is not None:
    lst= list(idjob.values())
    print(tabulate(lst, maxlen, asTsv))

