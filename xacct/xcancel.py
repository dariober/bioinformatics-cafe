#!/usr/bin/env python3

import sys
import subprocess
import csv
#import collections
#from io import StringIO
import argparse
import os
#import datetime
import re

STATES= ['PENDING', 'RUNNING', 'SUSPENDED']

parser = argparse.ArgumentParser(description= """
Cancel slurm jobs by convienient matching of job names.

EXAMPLES

xcancel.py -n myJob           # cancel jobs matching this name
xcancel.py -n '.*myJob$'      # cancel jobs with name ending in 'myJob'
xcancel.py -n myJob -t RUN PD # Appy to jobs in RUNNING or PENDING state

""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--name', '-n',
                   default= '',
                   help='''Cancel jobs with name matching this regular expression. MEMO: wrap pattern in
single quotes to avoid shell expansion.
                   ''')

parser.add_argument('--state', '-t',
                   nargs= '*',
                   default= STATES,
                   help='''Cancel jobs in these states. State can be abbreviated and is case insensitive.
Default %(default)s
                   ''')

parser.add_argument('--dryrun', '-d',
                   action= 'store_true',
                   help='''Only show the jobs that would be cancelled and the scancel command.
                   ''')

parser.add_argument('--verbose', '-V',
                   action= 'store_true',
                   help='''Verbose output for debugging only
                   ''')

parser.add_argument('scancel_args',
                   nargs= "*",
                   help='''Further args to scancel. Sperate them from the other arguments with `--`''')

parser.add_argument('--version', '-v', action='version', version='%(prog)s 0.1.0')

args= parser.parse_args()

# -----------------------------------------------------------------------------

states= set()
for x in args.state:
    x= x.upper()
    if x == 'PD':
        x= 'PENDING'
    found= False
    for choice in STATES:
        if choice.startswith(x):
            states.add(choice)
            found= True
            break
    if not found:
        sys.exit('Unrecognized state %s.' % x)
            
try:
    name_pattern= re.compile(args.name)
except:
    sys.stderr.write("Pattern '%s' does not compile to a valid regular expression\n" % args.name)
    sys.exit(1)

cmd= ['sacct', '--parsable2', '--format=JobID,JobName%1000,State']
if args.verbose:
    sys.stderr.write(' '.join(cmd) + '\n')
try:
    sacct= subprocess.check_output(cmd, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as exc:
    print(exc.output.decode().strip())
    print('Exit code %s' % exc.returncode)
    sys.exit(1)

reader = csv.DictReader(sacct.decode().split('\n'), delimiter='|')
sacct_out= []
for line in reader:
    sacct_out.append(line)

jobs= []
for job in sacct_out:
    if job['State'] not in states:
        continue
    m= re.fullmatch(name_pattern, job['JobName'])
    if m is not None and job not in jobs:
        jobs.append(job)

cmd= ['scancel'] + args.scancel_args + [x['JobID'] for x in jobs]

header= 'JobID\tJobName\tState'
if args.dryrun:
    print(header)
    for job in jobs:
        print('{JobID}\t{JobName}\t{State}'.format(**job))
    print(header)
    sys.stderr.write(' '.join(cmd) + '\n')
    sys.exit()
else:
    if args.verbose:
        sys.stderr.write(' '.join(cmd) + '\n')
    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        print(exc.output.decode().strip())
        print('Exit code %s' % exc.returncode)
        sys.exit(1)
