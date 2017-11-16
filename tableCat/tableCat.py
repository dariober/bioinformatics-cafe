#!/usr/bin/env python

import argparse
import sys
import glob
import os
import gzip
import re
import subprocess

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Concatenate files and add a column of file identifiers. The output is
    going to be in "database long format".
    
    Typical use case: You have a pipeline that generates file in tabular format,
    one file per input, all files the same format. You want to concatenate these
    tables while keeping the identity of the files and optionally add a header
    for the cat'd files.

EXAMPLE
    cd tableCat
    ./tableCat.py -i test/*.txt -H
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--input', '-i',
    required= True,
    nargs= '+',
    help='''Input files to concatenate. Use - to read the list of files from stdin.
''')

parser.add_argument('--sep', '-s',
    default= '\t',
    help='''Separator used to right-bind the column with file id. Default to tab.
''')

parser.add_argument('--keepdir', '-kd',
    action= 'store_true',
    help='''Keep the the directory name in the output. Default is to strip it.
''')

parser.add_argument('--regex', '-r',
    default= None,
    help='''an optional regex to strip from the file name to make the ID nicer.
''')

parser.add_argument('--idColName', '-id',
    default= 'file_id',
    help='''Name for the additional column with the file id. Default 'file_id'.
''')

parser.add_argument('--skip', '-S',
    type= int,
    default= 0,
    help='''Skip this many lines from *each* file.
''')

parser.add_argument('--firsthdr', '-H',
    action= 'store_true',
    help='''First non-skipped line from *first* file is header and it is printed.
    The first non-skipped line from following files is NOT printed.
''')


parser.add_argument('--version', action='version', version='%(prog)s 0.2.0')

# -----------------------------------------------------------------------------

def which(program):
    """Test if program is on PATH, same as unix which.
    From http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def globToList(glist):
    """Expand glob(s) to list of individual files. Duplicates removed!
    """
    if type(glist) != list:
        sys.exit('Invalid arg')
    globbed= []
    for g in glist:
        xglobbed= glob.glob(g)
        for x in xglobbed:
            globbed.append(x)
    globbed= sorted(list(set(globbed)))
    return globbed

def isGzip(x):
    """Test if file x is gzip by attempting to open it with gzip module
    """
    try:
        fin= gzip.open(x, 'rb')
        fin.readline()
        fin.close()
        return True
    except IOError:
        return False


def openGzip(x):
    """Try to open file x with either system gunzip or with
    python gzip module. Return None if file is not gzip
    """
    if not isGzip(x):
        print("Not gzip")
        return None

    if which('gunzip') is None:
        fin= gzip.open(x, 'rb')
        fin.readline()
        fin.seek(0)
        return fin
    else:
        p= subprocess.Popen(['gunzip', '-c', x], shell= False, stdout= subprocess.PIPE, stderr= subprocess.PIPE)

        #err= []
        #for line in p.stderr:
        #    err.append(line)
        #
        #if len(err) != 0:
        #    raise Exception('%s\n' %('\n'.join(err)))
        return p.stdout

def xopen(x):
    """Seemlessly open flat or gzip file.
    x: File name
    Return: Open file handle ready to iterate through
    """
    if not isGzip(x):
        fin= open(x)
    else:
        fin= openGzip(x)
    return fin

def parseID(filename, keepdir= False, regex= None):
    """Parse the filename to optionally remove dir and a regex
    """
    id= filename
    if not keepdir:
        id= os.path.basename(filename)
    if regex:
        id= re.sub(regex, '', id)
    return(id)

if __name__ == '__main__':

    # Prevent BrokenPipe error when tableCat is piped to e.g. head. 
    # see http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL) 

    args= parser.parse_args()
    if args.input == ['-']:
        files= []
        for x in sys.stdin:
            files.append(x.strip())
    else:
       files= globToList(args.input)

    if files == []:
        sys.exit('No file found to concatenate!')
    is_first_file= True
    for f in files:
        # sys.stderr.write('File: %s ...' %(f))
        fin= xopen(f)
        id= parseID(f, args.keepdir, args.regex)
        is_first_line= True
        n= 0
        skip= args.skip
        for line in fin:
            line2= line.rstrip('\n') + args.sep + id
            if skip > 0:
                skip-=1
                continue
            if args.firsthdr:
                if is_first_file and is_first_line:
                    print(line.rstrip('\n') + args.sep + args.idColName)
                    n+=1
                else:
                    if is_first_line:
                        pass
                    else:
                        print(line2)
                        n+=1
                is_first_line= False
            else:
                print(line2)
                n+=1
        is_first_file= False
        # sys.stderr.write(' %s lines printed\n' %(n))
    sys.exit()
    
