#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import tempfile
import re
import gzip
import shlex

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Execute QC of BS/oxBS libraries.

PIPELINE
    - Shorten reads oxbs_qc.trim_fastq()
    - Remove adapters trim_galore
    - Map reads w/ bismark
    - Convert sam to bam, sort & index
    - Clena names & Clip overallaping reads `bam clipOverlap`
    - Call methylation 
    - Plot coverage
    - 

REQUIRES
    - fastx_trimmer
    - cutadapt
    - trim_galore
    - Bismark
    - Bowtie2
    - samtools
    - cleanBamReadNames.py -b -i
    - bamUtil clipOverlap
    - mpileup2methylation.py
    
bismark --bowtie2 -o bismark_out /lustre/sblab/berald01/reference_data/genomes/bsseq_synthetic4/ fastq_trimmed/$trimmedfq &&
mv bismark_out/${trimmedfq}_bt2_bismark.sam bismark_out/${bname}.bs4.sam &&
sam2bam.py bismark_out/${bname}.bs4.sam &&
rm fastq_trimmed/$trimmedfq

    """, formatter_class= argparse.RawTextHelpFormatter)

# -----------------------------------------------------------------------------
input_args= parser.add_argument_group('Input options', '')

input_args.add_argument('--fastq', '-f',
                   required= True,
                   nargs= '+',
                   help='''One or two fastq files to QC. If two files are passed
they will be treated as paired-ends
                   ''')

input_args.add_argument('--refdir', '-r',
                   required= True,
                   help='''The path to the directory containing the reference
sequences. This directory is expected to contain
- The reference FASTA file <ref>.fa
- The tab delimited file of base modifications <ref>.txt
- The subdir `Bisulfite_Genome` prepared by `bismark_genome_preparation --bowtie2`
                   ''')

input_args.add_argument('--outdir', '-o',
                   required= False,
                   default= 'oxbs_qc_out',
                   help='''Directory where output will be sent. It will be created
if it doesn't exist.
                   ''')

input_args.add_argument('--skip_shorten', '-Ss',
                   action= 'store_true',
                   help='''Do not perform shortening of reads. Use this option
if your reads are 3 or more bases shorter than the shortest reference sequence.
                   ''')

input_args.add_argument('--skip_trim', '-St',
                   action= 'store_true',
                   help='''Do not perform trim reads by removing 3'-adapters and
low quality ends.
                   ''')

input_args.add_argument('--bin',
                   required= False,
                   default= '',
                   help='''Path to the executable files required by oxbs_qc. Default 
is '' (all executables are on $PATH)
                   ''')

class FastqKit:
    """Given input fastq files, get the file names that will be associated with the
    various steps.
    """
    def __init__(self):
        self.arg_fq = []     ## List of 1 or 2 fastq files exactly as passed on args.fastq
        self.is_paired= None ## 
        self.is_gzip= None
        self.fq = []         ## List of 1 or 2 with the name (no path) of fastq files
        self.fq_path= []     ## List of 1 or 2 with the path only of fastq files
        self.fq_ext= ''         ## Extension of fastq files. Extesions must be the same for both e.g. '.fq.gz' '.fastq.gz'
        self.bname= []          ## Basename for the file. I.e. fastq name w/o path and extensions .txt, .fq, .fastq, .gz
        self.fastx_trimmer_fq= [] ## Path & Names of fastq files after fastx_trimmer
        self.trimg_fq= []         ## Path & Names of the fastq files after trim_galore
        self.allowed_ext= ('.txt', '.fq', '.fastq') ## Allowed extensions for fq files. Ignoring possible .gz
    def __str__(self):
        xl= []
        for k in sorted(self.__dict__):
            v= self.__dict__[k]
            xl.append(str(k) + ': ' + str(v))
        return('\n' + '\n'.join(xl) + '\n')
    def __eq__(self, other): 
        """This is to check equality of objects
        """
        return self.__dict__ == other.__dict__        
    
    def dissect_fastqname(self):
        """Populate slots in FastqKit with characteristincs if FASTQ file
        """
        if self.arg_fq == []:
            return(self)
        if len(self.arg_fq) == 2:
            self.is_paired= True
        if len(self.arg_fq) == 1:
            self.is_paired= False
        fq_ext= []
        self.fq_path= []
        self.fq= []
        self.bname= []
        for fq in self.arg_fq:
            p, f= os.path.split(fq)
            self.fq_path.append(p)
            self.fq.append(f)
            name, ext= os.path.splitext(re.sub('\.gz$', '', f)) ## Extension ignoring .gz
            fq_ext.append(ext)
            self.bname.append(name)          ## Basename: Fastq name w/o path, '.gz', extension
        if self.fq[0].endswith('.gz'):
            self.is_gzip= True
        else:
            self.is_gzip= False
        ## Check consistency
        if len(self.arg_fq) > 2:
            raise Exception('\nToo many fastq files!\n')
        invalid_ext= [x for x in fq_ext if x not in self.allowed_ext]
        if invalid_ext != []:
            raise Exception('Invalid extension %s' %(fq_ext))
        if self.is_paired:
            ## Consistency for paired-end
            if len(set(self.fq)) != 2:
                raise Exception('\nFastq file names must have different names, even if they are in different directories: %s\n' %(self.fq))
            if len(set(fq_ext)) != 1:
                raise Exception('\nFastq files must have the same extension: %s' %(fq_ext))
            gz= 0
            for f in self.fq:
                if f.endswith('.gz'):
                    gz += 1                    
            if gz != 0 and gz != 2:
                print(self)
                raise Exception('\nFastq files must be *both* either gzipped compressed uncompressed: %s\n' %(self.fq))
        self.fq_ext= fq_ext[0] ## Extension has been validated. Assign it.        
        return(self)
           
    def fastqname_to_trim_galore_name(self, dont_gzip= False, gzip= True):
        """Return the name of the fastq files that would be produced by trim_galore
        dont_gzip & gzip:
            same options that will be passed to trim_galore.
        Return:
            List of len 1 or 2 with the names that trim_galore will produce. Path
            is never included as this depends on trim_galore -o 

        Combinations are
            - gz: Yes/No
            - Ext: .fq / .fastq / .txt
            - PE: Yes/No
            
            SE
            mjb042_oxBS_R1.fastq.gz: mjb042_oxBS_R1_trimmed.fq.gz
            mjb042_oxBS_R1.fq.gz:    mjb042_oxBS_R1_trimmed.fq.gz
            mjb042_oxBS_R1.txt.gz:   mjb042_oxBS_R1.txt.gz_trimmed.fq.gz
            PE
            mjb042_oxBS_R1.fq.gz     mjb042_oxBS_R1_val_1.fq.gz
            mjb042_oxBS_R2.fq.gz     mjb042_oxBS_R2_val_2.fq.gz
            
            mjb042_oxBS_R1.fastq.gz     mjb042_oxBS_R1_val_1.fq.gz
            mjb042_oxBS_R2.fastq.gz     mjb042_oxBS_R2_val_2.fq.gz
            
            mjb042_oxBS_R1.txt.gz     mjb042_oxBS_R1.txt.gz_val_1.fq.gz
            mjb042_oxBS_R2.txt.gz     mjb042_oxBS_R2.txt.gz_val_2.fq.gz
            
            Output is always going to *.gz since we run trim_galore --gzip
        """
        if self.fq == [] or not self.fq:
            return(self)
        ext= self.fq_ext
        if not self.is_paired:
            if ext in ('.fastq', '.fq'):
                trim_g= [self.bname[0] + '.fq']
            else:
                trim_g= [self.bname[0] + self.fq_ext]
        else:
            if ext in ('.fastq', '.fq'):
                trim_g= [self.bname[0] + '_val_1.fq', self.bname[1] + '_val_2.fq']
            else:
                trim_g= [self.fq[0] + '_val_1.fq', self.fq[1] + '_val_2.fq']
        if gzip:
            self.trimg_fq= [x + '.gz' for x in trim_g]
        else:
            self.trimg_fq= trim_g
        return(self)

def trim_fastq(fq_in, mlen, fq_out):
    """Replace fastx_trimmer: Trim reads from the beginning to lenght given in
    mlen.
    fq_in, fq_out:
        Name of fastq files for input and output.
    mlen:
        Int
    Return:
        Possibly: FastqKit object.
        Side effect: Fastq file trimmed and gzipped.
    """
    if fq_in.endswith('.gz'):
        fin= gzip.open(fq_in)
    else:
        fin= open(fq_in)
    fout= gzip.open(fq_out, 'wb')
    n= -1
    for line in fin:
        if n % 2 == 0:
            line= line.rstrip()[0:mlen] + '\n'
        fout.write(line)
        n += 1
    fin.close()
    fout.close()   
    
def getReferenceSize(refdir):
    """Get the length shortest sequence from the (control) fasta file of control sequences
    refdir:
        Directory as it would be passed to bismark
    Return:
        Int (length of shortes seq)
    """
    fastalist= os.listdir(refdir)
    fasta= [x for x in fastalist if x.endswith('.fa')]
    if len(fasta) != 1:
        raise Exception('One and only one fasta file must be found in %s' %(refdir))
    refseq= os.path.join(refdir, fasta[0])
    try:
        fai= open(refseq + '.fai')
    except:
        cmd= 'samtools faidx %s' %(refseq)
        p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
        stdout, stderr= p.communicate()
        if stderr != '':
            print(cmd)
            print(stderr)
            raise Exception('Failed to index %s' %refseq)
        fai= open(refseq + '.fai')
    idx= fai.readlines()
    idx= [x.strip().split('\t') for x in idx]
    idx= [int(x[1]) for x in idx]
    shortest= min(idx)
    return(shortest)

def get_fastq_encoding(fastq):
    """Read fastq file to determine the quality encoding. Fastq file can be either
    plain text or gzipped.
   
    Returns:
        String ``Sanger`` or ``Illumina 1.5+`` or ``Undetermined``.
        There is no distinction within Illumina variations and no distictions between Sanger and `Illumina 1.8+ Phred+33`.
        If only ambiguos codes are found (e.g. 'BCDEFGHI') the string ``Undertermined`` is returned
    """
    if fastq.endswith('.gz'):
        fin= gzip.open(fastq)
    else:
        fin= open(fastq)
    while True:
        fin.readline()
        fin.readline()
        fin.readline()
        qual= fin.readline().strip()
        if qual == '':
            fin.close()
            return('Undetermined')
            ## sys.exit(inspect.stack()[0][3] + ': Encoding of %s could not be determined after having read the whole file' %(fastq))
        encoding= [ord(x) for x in qual]
        for x in encoding:
            if x < 59:
                fin.close()
                return('Sanger')
            if x > 74:
                fin.close()
                return('Illumina 1.5+')

def get_wdir(wdir):
    """Set a working directory.
    Possibly TODO: use tempfile.mkdir if wdir is None 
    """
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    return(wdir)


def task_trim_fastq(inFastqKit, refsize, outdir):
    """Execute trim_fastq on the FastqKit object.
    refsize:
        int for length of reference. Trimming will be refsize - 3
    outdir:
        Output directory.
    return:
        FastqKit object with trimmed seq
    """
    outdir= get_wdir(outdir)
    fq_short= []
    for i in range(len(inFastqKit.arg_fq)):
        fqx= inFastqKit.bname[i] + '.x' + inFastqKit.fq_ext + '.gz'
        fqx= os.path.join(outdir, fqx)
        trim_fastq(inFastqKit.arg_fq[i], refsize - 3, fqx)
        fq_short.append(fqx)
    oxbs_short= FastqKit()
    oxbs_short.arg_fq= fq_short
    oxbs_short.dissect_fastqname()
    return(oxbs_short)

def task_trim_galore(inFastqKit, opts= '', path= ''):
    """Execute trim_galore on the input FastqKit.
    inFastqKit:
        FastqKit object
    opts:
        String of options passed to trim_galore
    path:
        path to trim_galore
    return:
        Tuple with 1) FastqKit object with the output of trim_galore has `arg_fq`
        2) the command line passed to the system call as *list*.
    """
    outFastqKit= FastqKit()
    inFastqKit.dissect_fastqname()
    opts= shlex.split(opts)
    if '-o' in opts:
        i= opts.index('-o')
        outdir= opts[i+1]
    elif '--output_dir' in opts:
        i= opts.index('--output_dir')
        outdir= opts[i+1]
    else:
        outdir= ''
    if '--paired' in opts and not inFastqKit.is_paired:
        print(inFastqKit.is_paired)
        raise('--paired opt passed to trim_galore but fastq files are not marked as paired!')
    elif '--paired' not in opts and inFastqKit.is_paired:
        opts.append('--paired')
    else:
        pass
    trim_galore= os.path.join(path, 'trim_galore')

    cmd= [trim_galore] + opts + inFastqKit.arg_fq
    print(' '.join(cmd))
    p= subprocess.Popen(' '.join(cmd), shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        raise Exception('Error produced by trim_galore')
    inFastqKit.fastqname_to_trim_galore_name()
    outFastqKit.arg_fq= inFastqKit.trimg_fq
    outFastqKit.dissect_fastqname()
    outFastqKit.fq_path= [outdir] * len(inFastqKit.trimg_fq)
    print('Trimmed fastq(s) path: %s' %(', '.join(set(outFastqKit.fq_path))))
    print('Trimmed fastq(s) name: %s' %(', '.join(outFastqKit.arg_fq)))
    return((outFastqKit, cmd))

def task_bismark_aln(inFastqKit, genome_folder, opts= '', path= ''):
    """Execute bismark on the input FastqKit.
    inFastqKit:
        FastqKit object
    opts:
        String of options passed to bismark
    path:
        path to bismark
    return:
        Tuple with 1) Name of output bam sam file
        2) the command line passed to the system call as *list*.
    """
    opts= shlex.split(opts)
    if '-1' in opts or '-2' in opts:
        print('bismark options -1 and -2 are not allowed. Pass fastq files via `oxbs_qc.py --fastq`')
        raise Exception('Not allowed option passed to bismark')
    if '-o' in opts:
        i= opts.index('-o')
        outdir= opts[i+1]
    elif '--output_dir' in opts:
        i= opts.index('--output_dir')
        outdir= opts[i+1]
    else:
        outdir= ''
    if inFastqKit.is_paired:
        fq1= os.path.join(inFastqKit.fq_path[0], inFastqKit.fq[0])
        fq2= os.path.join(inFastqKit.fq_path[1], inFastqKit.fq[1])
        fq= '-1 %s -2 %s' %(fq1, fq2)
    else:
        fq= os.path.join(inFastqKit.fq_path[0], inFastqKit.fq[0])
    if '--bowtie2' not in opts:
        bwt2= '--bowtie2'
    else:
        bwt2= ''

    bismark= os.path.join(path, 'bismark')

    cmd= [bismark] + opts + fq
    print(' '.join(cmd))
    p= subprocess.Popen(' '.join(cmd), shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        raise Exception('Error produced by trim_galore')    
    
def main():
    args = parser.parse_args()
    wdir= get_wdir(args.outdir)
    encoding= get_fastq_encoding(args.fastq[0])
    refsize= getReferenceSize(args.refdir)
    fqkit= FastqKit()
    fqkit.arg_fq= args.fastq
    fqkit.dissect_fastqname()
    if not args.skip_shorten:
        print('Shortening reads: %s' %(', '.join(fqkit.arg_fq)))
        fqkit= task_trim_fastq(fqkit, refsize, wdir)
        print('Shortened reads in: %s' %(', '.join(fqkit.arg_fq)))
    ## trim galore
    if not args.skip_trim:
        fqkit= task_trim_galore(fqkit, opts= '-o %(wdir)s' %{'wdir': wdir}, path= args.bin)
    ## Bismark alignment:
    
    sys.exit()
        
if __name__ == '__main__':
    main()
    
