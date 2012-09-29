#!/usr/bin/env python

import argparse
import sys
import subprocess
import md5
import sys
import sblab
import os
import shutil
import shlex

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Execute fastqc and computes md5sum of the input files. Md5sums are are written
    to <input.fq>_fastqc/fastqc_data.txt as last line of the "Basic Statistics" module.
    
EXAMPLE

TODO:

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fastq', '-i',
                   required= True,
                   nargs= '+',
                   help='''List of input files to pass to fastqc and to compute
md5sum.''')

parser.add_argument('--fastqc', '-f',
                   required= False,
                   default= '',
                   type= str,
                   help='''String of arguments to pass as is to fastqc. E.g. ' --noextract -t 8'.
(The list of input files not included.). Default is ''.
NB: Use quotes and start with a blank space: E.g.
OK:    ' --noextract'
WRONG: '--noextract' 
                   ''')

parser.add_argument('--fastqc_path', '-p',
                   type= str,
                   required= False,
                   default= '',
                   help='''Path to fastqc. Default is '' meaning fastqc is on the
path.
                   ''')

args = parser.parse_args()

## ----------------------------------------------------------------------------

## {{{ http://code.activestate.com/recipes/266486/ (r1)
def sumfile(fobj):
    '''Returns an md5 hash for an object with read() method.'''
    m = md5.new()
    while True:
        d = fobj.read(8096)
        if not d:
            break
        m.update(d)
    return m.hexdigest()


def md5sum(fname):
    '''Returns an md5 hash for file fname, or stdin if fname is "-".'''
    if fname == '-':
        ret = sumfile(sys.stdin)
    else:
        try:
            f = file(fname, 'rb')
        except:
            return 'Failed to open file'
        ret = sumfile(f)
        f.close()
    return ret

def add_md5_fastqc(fastqc_data_file, md5):
    """Add md5sum to fastqc_data.txt as last line of Basic Statistics module.
    Return a list where each line is a line of the input.
    """
    fastqc_data= open(fastqc_data_file).readlines()    
    n= 0
    for line in fastqc_data:
        if line.startswith('>>END_MODULE'):
            break
        else:
            n += 1
    fastqc_data.insert(n, 'md5sum\t%s\n' %(md5))
    return(fastqc_data)
# -----------------------------------------------------------------------------

# Parse fastqc cmd option to get output dir:
fastqccmd= shlex.split(args.fastqc)
if '-o' in fastqccmd and '--outdir' in fastqccmd:
    sys.exit('Invalid commands passed to fastqc: "%s".' %(args.fastqc))
elif '-o' in fastqccmd:
    outd= os.path.abspath(fastqccmd[fastqccmd.index('-o') + 1])
elif '--outdir' in fastqccmd:
    outd= os.path.abspath(fastqccmd[fastqccmd.index('--outdir') + 1])
else:
    outd= None

## Remove possible duplicates
fastq= [os.path.abspath(x) for x in args.fastq]
fastq= sorted(set(fastq))

n= 0
for f in fastq:
    try:
        fh= open(f)
        fh.close()
    except:
        print('Skipping file %s: File not found' %(f))
        del fastq[n]
    n += 1
if fastq == []:
    sys.exit('\nNo file found- Nothing to be done!\n')
    
fastqccmd= os.path.join(args.fastqc_path, 'fastqc') + args.fastqc + ' ' + ' '.join(fastq)
print(fastqccmd)
p= subprocess.Popen(fastqccmd, shell= True)
p.wait()

for fq in fastq:
    "get md5sums"
    if outd is None:
        """fqdir is where the output files from fastqc should be found. This is
        the same dir where inputs are if --outdir was not specified"""
        fqdir= os.path.split(fq)[0]
        if fqdir == '':
            fqdir= './'
    else:
        "Output fastqc files are where --outdir was set"
        fqdir= outd
    print('Computing md5sum for: %s' %(fq))
    md5= md5sum(fq)
    fastqc_dir= os.path.join(fqdir, os.path.split(sblab.get_fastqc_dir(fq))[1])
    cmd= 'unzip -q -o -d %s %s.zip' %(fqdir, fastqc_dir)
    print(cmd)
    p= subprocess.Popen(cmd, shell= True)
    p.wait()
    fastqc_data_file= os.path.join(fastqc_dir, 'fastqc_data.txt')
    fastqc_data= add_md5_fastqc(fastqc_data_file, md5)
    fastqc_data_out= open(os.path.join(fastqc_dir, 'fastqc_data.txt'), 'w')
    for line in fastqc_data:
        "Replace original file"
        fastqc_data_out.write(line)
    fastqc_data_out.close()
    cmd= 'zip -q %s.zip %s' %(fastqc_dir, fastqc_dir)
    print(cmd)
    p= subprocess.Popen(cmd, shell= True)
    p.wait()
    if ' --noextract ' in args.fastqc:
        shutil.rmtree(fastqc_dir)
sys.exit()