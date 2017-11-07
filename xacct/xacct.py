#!/usr/bin/env python3

import sys
import subprocess
import csv
import collections
from io import StringIO
import argparse
import os
import datetime

parser = argparse.ArgumentParser(description= """
Execute slurm sacct and return a friendly tabular output. 
Values in MaxRSS and ReqMem are in MB units.

EXAMPLES

xacct.py
xacct.py -- -S 2016-12-01  # Show jobs starting from YYYY-MM-DD
xacct.py -d 1              # Jobs from yesterday onwards

# Sort by memory usage
xacct.py -tsv | tail -n+2 | sort -t'   ' -k4,4n
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--days', '-d',
                   type= float,
                   default= 0,
                   help='''Show jobs from this many days ago. Deafult %(default)s''')

parser.add_argument('--fromId', '-id',
                   type= int,
                   default= None,
                   help='''Show jobs whose ID is equal or greater than this. Deafult %(default)s''')

parser.add_argument('--tsv', '-tsv',
                   action= "store_true",
                   help='''Print columns separated by TAB (better for further processing) 
instead of tabulating them using spaces (better for eyeballing). This option automatically sets also
--no-color''')

parser.add_argument('--no-color', '-nc',
                   action= "store_true",
                   dest= 'no_color',
                   help='''Do not add color to the output strings. Use this option if you need
to parse the output and the color codes strings get on your way.
                   ''')

parser.add_argument('--verbose', '-V',
                   action= "store_true",
                   help='''Verbose for debugging. Print to stderr the sacct command.''')

parser.add_argument('sacct_args',
                   nargs= "*",
                   help='''Further args to sacct, e.g. `-S 2017-09-10`. Do not include `--format`.
Use `--` to add slurm args.''')

parser.add_argument('--version', '-v', action='version', version='%(prog)s 0.1.0')

args= parser.parse_args()

# -----------------------------------------------------------------------------

def normalizeMem(x):
    """We use binary multiplier instead of decimal to convert kilo and giga to
    mega. I.e.  1024K = 1M. Compare to help for --noconvert option: `Don't
    convert units from their original type (e.g. 2048M won't be converted to
    2G).`. 

    normalizeMem('2G') -> 2048
    """
    if x.strip() == '':
        return '';  
    x= x.strip('n').strip('c') # ReqMem appends 'n' or 'c', see `man sacct` 
    mem= -1
    if x.endswith('K'):
        mem= float(x.strip('K')) * 2**10
    elif x.endswith('M'):
        mem= float(x.strip('M')) * 2**20
    elif x.endswith('G'):
        mem= float(x.strip('G')) * 2**30
    else:
        mem= float(x)
    return round(mem/(2**20))

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

def tabulate(line, maxlen, asTsv, no_color):
    if asTsv:
        return '\t'.join([str(x) for x in line])
    if not no_color:
        line= colorize(line)
    row_fmt= "{:<8}{:<%s}{:<16}{:<8}{:<8}{:<10}{:<12}{:<30}" % (maxlen + 2)
    return row_fmt.format(*line)

def colorize(lst):
    """Add color codes to the printable list
    For codes see https://misc.flogisoft.com/bash/tip_colors_and_formatting
    """
    xcol= []
    for x in lst:
        if x == 'FAILED':
            x= '\033[31mFAILED\033[0m'
        elif x == 'COMPLETED':
            x= '\033[32mCOMPLETED\033[0m'
        elif x == 'RUNNING':
            x= '\033[94mRUNNING\033[0m'
        xcol.append(x)
    return xcol

# -----------------------------------------------------------------------------

starttime= []
if args.days > 0:
    d= datetime.datetime.today() - datetime.timedelta(days= args.days)
    starttime= ['--starttime', d.date().isoformat()]

cmd= ['sacct', '--parsable2', '--format=JobID,JobName%50,NodeList,MaxRSS,ReqMem,AllocCPUS,Elapsed,State'] + starttime + args.sacct_args
if args.verbose:
    sys.stderr.write(' '.join(cmd) + '\n')
sacct= subprocess.check_output(cmd)

reader = csv.DictReader(sacct.decode().split('\n'), delimiter='|')

maxlen= maxNameLen(sacct)
if args.tsv:
    no_color= True
else:
    no_color= args.no_color

print(tabulate(reader.fieldnames, maxlen, args.tsv, no_color))

try:
    idjob= None
    for line in reader:
        if args.fromId is not None and args.fromId >= int(line['JobID'].replace('.batch', '')):
            continue
        line['MaxRSS']= normalizeMem(line['MaxRSS'])
        line['ReqMem']= normalizeMem(line['ReqMem'])
        if '.batch' not in line['JobID']:
                if idjob is not None:
                    # This line is not a batch job, so the previous line
                    # can be printed as it doesn't have an associated batch job.
                    lst= list(idjob.values())
                    print(tabulate(lst, maxlen, args.tsv, no_color))
                idjob= line
        else:
            # This is a batch job. So associate to it the job ID line    
            if idjob['JobID'] == line['JobID'].replace('.batch', ''):
                lst= list(fillInJob(idjob, line).values())
                print(tabulate(lst, maxlen, args.tsv, no_color))
            else:
                print(idjob)
                print(line)
                raise Exception("Cannot find job id")
            idjob= None
    if idjob is not None: # and args.fromId > 0 and args.fromId >= int(line['JobID']):
        lst= list(idjob.values())
        print(tabulate(lst, maxlen, args.tsv, no_color))
    print(tabulate(reader.fieldnames, maxlen, args.tsv, no_color))

except (BrokenPipeError, IOError):
    # This is to avoid stack trace in e.g. `sacct_easy.py | head`
    # See also https://stackoverflow.com/questions/26692284/brokenpipeerror-in-python
    pass

