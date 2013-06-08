import sys
import os
import subprocess
import tempfile
import re
import gzip
import shlex
import shutil


class TrimGaloreException(Exception):
    """Base class for exceptions raised by trim_galore."""
    pass

class TrimFastqException(Exception):
    """Base class for exceptions raised by trim_galore."""
    pass

class BismarkException(Exception):
    pass

class GetWdirException(Exception):
    pass

class Sam2BamException(Exception):
    pass


def get_wdir(wdir):
    """Set a working directory.
    Possibly TODO: use tempfile.mkdir if wdir is None 
    """
    if wdir == '' or type(wdir) != str:
        raise GetWdirException('wdir must be a non-empty str. Got "%s"' %(wdir))
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    return(wdir)

def task_trim_galore(inFastq, outFastq, tmpdir= 'trim_galore_wdir', opts= '', path= ''):
    """Execute trim_galore on the input FastqKit.
    inFastq:
        List of length one or two. If length 2 it is taken as paired-end
    tmpdir:
        A working directory where to dump output files. Inside this dir python will
        make a tempdirectory and the ouput of trim_galore put there
    outFastq:
        List of length 1 or 2 for the name of the output fastq file(s). DO NOT include
        path (it will be stripped anyway). 
    opts:
        String of options passed to trim_galore
    path:
        path to trim_galore
    echo:
        Only print the command that would be executed.
    return:
        Dictionary with:
        {'cmd': 'command-executed-by-subprocess',
         'fastq': ['fastq1', 'fastq2'],
         'report': 'trim_galore-report'}
    """
    if type(inFastq) != list:
        raise TrimGaloreException('Input fastq files must be a *list*. Got %s' %(type(inFastq)))
    if type(outFastq) != list:
        raise TrimGaloreException('Output fastq files must be a *list*. Got %s' %(type(outFastq)))
    if len(inFastq) < 1 or len(inFastq) > 2:
        raise TrimGaloreException('Expected a list of 1 or 2 fastq files. Got %s.' %(len(inFastq)))
    for x in inFastq:
        if not os.path.exists(x):
            raise TrimGaloreException('File %s not found' %(x))
    for i in range(0, len(outFastq)):
        fq= outFastq[i]
        p,f= os.path.split(fq)
        if p != '':
            raise TrimGaloreException('Output fastq file(s) must not have directory path. Got "%s"' %(fq))
    if len(inFastq) != len(outFastq):
        raise TrimGaloreException('%s fastq files in input but %s found in output arg' %(len(inFastq), len(outFastq)))
    if len(set(inFastq)) != len(inFastq):
        raise TrimGaloreException('Invalid input (same file given twice?): %s' %(inFastq))
    if len(set(outFastq)) != len(outFastq):
        raise TrimGaloreException('Invalid input (same file given twice?): %s' %(outFastq))
    ## Set working dir:
    tg_outdir= tempfile.mkdtemp(prefix= 'tmp_trim_galore_', dir= get_wdir(tmpdir))
    opts= shlex.split(opts)
    
    ## This otuput dir will NOT be passed to trim_galore. Instead the output files
    ## from trim_galore will be moved there afterwards.
    if '-o' in opts:
        i= opts.index('-o')
        final_outdir= opts[i+1]
        del opts[i:(i+2)]
    elif '--output_dir' in opts:
        i= opts.index('--output_dir')
        final_outdir= opts[i+1]
        del opts[i:(i+2)]
    else:
        final_outdir= '.'
    
    ## Handling gzip
    ## If gzip was given to trim_galore add *.gz to output name.
    ## If --dont_gzip id given remove *.gz from output names.
    ## -----------------
    if '--gzip' in opts:
        outFastq= [x.rstrip('.gz') + '.gz' for x in outFastq]
    if '--dont_gzip' in opts:
        outFastq= [x.rstrip('.gz') for x in outFastq]
    
    get_wdir(final_outdir)
    outFastqFull= [os.path.join(final_outdir, x) for x in outFastq]

    if '--paired' not in opts and len(inFastq) > 1:
        opts.append('--paired')
    else:
        pass
    trim_galore= os.path.join(path, 'trim_galore')

    cmd= [trim_galore] + opts + ['--output_dir'] + [tg_outdir] + inFastq
    print(' '.join(cmd))
    p= subprocess.Popen(' '.join(cmd), shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        raise TrimGaloreException('Error produced by trim_galore')
    ## Now get the names of the output fastq file(s) and rename them to outFastq
    ## -------------------------------------------------------------------------
    tg_output_files= sorted(os.listdir(tg_outdir))
    fqz= [x for x in tg_output_files if x.endswith('.fq.gz')]
    fq= [x for x in tg_output_files if x.endswith('.fq')]
    ## Make sure you got the right number of files.
    if len(fqz) > 2 or len(fq) > 2:
        raise TrimGaloreException('Too many fastq files found in output: %s, %s' %(fqz, fq))
    if len(fqz) > 0 and len(fq) > 0:
        raise TrimGaloreException('Found both gzipped and unzipped fastq files: %s, %s' %(fqz, fq))
    if len(fqz) != len(inFastq) and len(fq) != len(inFastq):
        if len(fqz) == 0:
            nout= len(fq)
        else:
            nout= len(fqz)
        raise TrimGaloreException('Inconsistent number of output fastq files found: %s input, %s output' %(len(inFastq), len(nout)))
    tg_fq= [os.path.join(tg_outdir, x) for x in fqz + fq]
    for old, new in zip(tg_fq, outFastqFull):
        os.rename(old, new)

    ## Get report:
    ## -----------
    tg_report=       [x for x in tg_output_files if x.endswith('_trimming_report.txt')]
    tg_report_tmp=   [os.path.join(tg_outdir, x) for x in tg_report] 
    tg_report_final= [os.path.join(final_outdir, x) for x in tg_report]
    if len(tg_report) != len(inFastq):
        raise TrimGaloreException('%s txt files found in trim_galore output dir. Expected %s (one per input).' %(len(tg_report), len(inFastq)))
    for src,dest in zip(tg_report_tmp, tg_report_final):
        shutil.move(src, dest)

    ## Exit.
    ## -----
    if len(os.listdir(tg_outdir)) != 0:
        raise TrimGaloreException('Temporary output dir from trim_galore should be empty now! Found: %s' %(sorted(os.listdir(tg_outdir))))
    shutil.rmtree(tg_outdir)
    return({'fastq': outFastqFull, 'cmd': cmd, 'report': tg_report_final})

def trim_fastq(fq_in, mlen, fq_out):
    """Replace fastx_trimmer: Trim reads from the beginning to lenght given in
    mlen.
    fq_in, fq_out:
        Name of fastq files for input and output.
    mlen:
        Int
    Return:
        Side effect: Fastq file trimmed and gzip'd. NB: .gz extension always present in output file!.
    """
    fq_out= fq_out.rstrip('.gz') + '.gz'
    if not os.path.exists(os.path.split(fq_out)[0]) and os.path.split(fq_out)[0] != '':
        raise TrimFastqException('Directory to output file %s does not exist' %(os.path.split(fq_out)[0]))
    if mlen <= 0:
        raise TrimFastqException('Invalid lenght to trim: %s' %(mlen))
    if fq_in == fq_out:
        raise TrimFastqException('Input and output files are the same!')
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
    return(True)

def task_bismark_aln(inFastq, outSam, ref, tmpdir= 'bismark_wdir', opts= '', path= '', echo= False):
    """Execute bismark alignment.
    inFastq:
        List of len 1 or 2 of fastq files. If 2 they are taken as paired.
    outSam:
        Output name for output sam file. Do not inlcude path (it will be stripped).
        Must end in *.sam
    ref:
        Path to reference.
    tmpdir:
        Working dir for bismark. It is NOT where you find the output. Output is given
        by `bismark --output_dir`
    opt:
        String of options to be passed to bismark
    path:
        Path to bismark
    """
    if type(inFastq) != list:
        raise BismarkException('Input fastq files must be a *list*. Got %s' %(type(inFastq)))
    if len(inFastq) < 1 or len(inFastq) > 2:
        raise BismarkException('Expected a list of 1 or 2 fastq files. Got %s.' %(len(inFastq)))
    if len(set(inFastq)) != len(inFastq):
        raise BismarkException('Input fastq files are the same! Got %s.' %(inFastq))
    for x in inFastq:
        if not os.path.exists(x):
            raise BismarkException('File %s not found' %(x))
    if not os.path.exists(ref):
        raise BismarkException('Path to reference %s not found' %(ref))
    if not outSam.endswith('.sam'):
        raise BismarkException('Output sam file must have extension .sam. Got %s' %(outSam))
    if os.path.split(outSam)[0] != '':
        raise BismarkException('Output sam name must not inlcude directory (use opts= "bismark --output_dir ..."). Got %s' %(outSam))
    ## Set working dir:
    bis_outdir= tempfile.mkdtemp(prefix= 'tmp_bismark_wdir_', dir= get_wdir(tmpdir))
    opts= shlex.split(opts)
    if '--bowtie2' not in opts:
        opts.append('--bowtie2')

    bismark= os.path.join(path, 'bismark')

    ## This otuput dir will NOT be passed to bismark. Instead the output files
    ## from bismark will be moved there afterwards.
    if '-o' in opts:
        i= opts.index('-o')
        final_outdir= opts[i+1]
        del opts[i:(i+2)]
    elif '--output_dir' in opts:
        i= opts.index('--output_dir')
        final_outdir= opts[i+1]
        del opts[i:(i+2)]
    else:
        final_outdir= '.'

    ## Single or paired end
    if len(inFastq) == 1:
        fq= [inFastq[0]]
    else:
        fq= ['-1', inFastq[0], '-2', inFastq[1]]
    cmd= [bismark] + ['-o'] + [bis_outdir] + opts + [ref] + fq
    cmd= ' '.join(cmd)
    if not echo:
        p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
        stdout, stderr= p.communicate()
        if p.returncode != 0:
            print(stderr)
            print('Exit code %s' %(p.returncode))
            raise BismarkException('Error produced by bismark')
    ## Move output. Sam will be renamed. Everything else stays the same
    ## ----------------------------------------------------------------
    output_files= sorted(os.listdir(bis_outdir))
    output_files_path= [os.path.join(bis_outdir, x) for x in output_files]
    sam_output= [x for x in output_files if x.endswith('.sam')]
    if len(sam_output) != 1:
        raise BismarkException('Expected 1 sam output. Found %s (%s)' %(len(sam_output), sam_output))
    sam_output= sam_output[0]
    sam= os.path.join(final_outdir, outSam)
    shutil.move(os.path.join(bis_outdir, sam_output), sam)
    all_out= [sam]
    for x in output_files_path:
        if x.endswith('.sam'):
            continue
        dst= os.path.join(final_outdir, os.path.split(x)[1])
        shutil.move(x, dst)
        all_out.append(dst)
    os.removedirs(bis_outdir)
    return({'cmd': cmd, 'sam': sam, 'all_out': all_out})

def cleanSamReadNames(inSam):
    """Removes trailing /1 /2 from read names to have pairs with the same name
    Output will replace input.
    """
    if not os.path.exists(inSam):
        raise Sam2BamException('Could not find sam file "%s"' %(inSam))
    if not inSam.endswith('.sam'):
        raise Sam2BamException('Invalid extension in %s. Expected .sam' %(inSam))
    samfile= open(inSam)
    tmpout= tempfile.NamedTemporaryFile(prefix= 'tmp_sam', suffix= '.sam', delete= False)
    n= 0
    header= True
    for read in samfile:
        if read.startswith('@'):
            tmpout.write(read)
            continue
        if n % 4 == 0:
            samline= read.split('\t')
            readName= samline[0]
            if readName[-2:] in ('/1', '/2'):
                readName= readName[0:-2]
            samline[0]= readName
            tmpout.write('\t'.join(read))
        else:
            tmpout.write(read)
        n += 1
    samfile.close()
    tmpout.close()
    shutil.move(tmpout.name, inSam)
    return(True)
    
def task_sam2bam(inSam, rmSam= True, path= ''):
    """Convert sam to bam, sort and index
    inSam:
        sam file to convert
    rmSam:
        Remove original SAM?
    path:
        Path to samtools
    return:
        Dictionary with list of cammands executed ('cmd_list'), output bam file ('bam')
    """
    if not os.path.exists(inSam):
        raise Sam2BamException('Could not find sam file "%s"' %(inSam))
    if not inSam.endswith('.sam'):
        raise Sam2BamException('Invalid extension in %s. Expected .sam' %(inSam))
    outBamUnsorted= re.sub('\.sam$', '.unsorted.bam', inSam)
    outBamPrefix= re.sub('\.sam$', '', inSam)
    outBam= outBamPrefix + '.bam'
    samtools= os.path.join(path, 'samtools')
    cmds= []
    cmds.append(' '.join([samtools, 'view', '-h', '-S', '-b', '-o', outBamUnsorted, inSam]))
    cmds.append(' '.join([samtools, 'sort',  outBamUnsorted, outBamPrefix]))
    cmds.append(' '.join([samtools, 'index', outBam]))
    for cmd in cmds:
        p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
        stdout, stderr= p.communicate()
        if p.returncode != 0:
            print(stderr)
            print('Exit code %s' %(p.returncode))
            print(cmd)
            raise Sam2BamException('Error produced while converting sam to bam')
    os.remove(outBamUnsorted)
    if rmSam:
        os.remove(inSam)
    return({'cmd_list': cmds, 'bam': outBam})

def task_clipOverlap(inBam, clippedBam):
    """Excute clipOverlap
    inBam:
        Input bam file to clip. Must be sorted.
        See http://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap
    clippedBam:
        Output clipped bam
    
    bam clipOverlap --in <inputFile> --out <outputFile>    
    """
    return()
#    cleanBamReadNames.py -b -i ${bismarksam%.sam}.bam | bam clipOverlap  --stats --in -.bam --out bismark_out/$bname.mm9.bam &&
#    samtools index bismark_out/$bname.mm9.bam