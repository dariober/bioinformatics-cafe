#!/usr/local/bin/python

import os
import subprocess
import sys
import shutil
import argparse
import getpass

parser = argparse.ArgumentParser(description= """

DESCRIPTION:

    Extract sequences from a fasta file corresponding to the input bed file
    (or any file with genomic regions) and run meme for motif enrichment analysis.
    Input bed can be pre-processed to filter and select regions.

    This program is an interface to the programs (in order of execution):
    
        awk  :         Filter input bed file
        sort :         Sort the filtered bed file
        head :         Use top n regions
        extend_bed :   Extends bed coordinates (internal function)
        FastaFromBed : Extract sequences corresponding to selected regions
        RepeatMasker : Pass sequence to it
        meme suite :   Find motifs
        tar/scp :      Tar and ship results to remote host
    
    Make sure the input file is in bed format (or convert it to bed using the
    --awk arg.)

    The output is sent to dir 'meme' created (if it doesn't exist) in the
    current dir, subdir input. E.g. if input bed is myfile.bed the output will
    be in ./meme/myfile This dir will be deleted at start of the run if it
    already exists.

    The order of the options below reflects the order of the input/output
    processing.
    
EXAMPLES:
   
    master_chipseq_motif.py 
        -i macs/peaks/el001_6.macs_combined.bed  # Input bed (without header)
        --awk '$14 < 5'                          # Include rows where 14th col is < 5
        --sort '-k 12,12nr -k 5,5nr'             # Sort by cols 12th and 5th (desc.) and...
        --npeaks 100                             # ...get top 100 rows
        --extendby 100                           # Extend bed coordinates by 100 bp left and right
        --scp '143.65.172.155:~/Tritume/'        # Send results to this remote dir
    
MEMO:
    The combined output from macs has no header and the columns are:

    *.macs_combined.bed:
    ------------
    1.  chr
    2.  summit_start
    3.  summit_end
    4.  macs_id
    5.  pileup_at_summit
    6.  strand
    7.  peak_start
    8.  peak_end
    9.  peak_length
    10. summit_dist_from_start
    11  tags_in_region
    12. log_pvalue  = -10*log10(pvalue)
    13. fold_enrich
    14. fdr_percent
    15. file_basename
    16. project_id
           

REQUIREMENTS:

    - bash tools awk and sort
    - bedtools' FastaFromBed
    - RepeatMasker suite (required if you want to mask repeats)
    - Meme suite

TODO:
    - Add options for the auxilary programs of the meme suite (tomtom, dreme etc.).
      Currently only the meme and fasta-subsample programs take options
    - Parameter strings starting with '-' are mis-interpreted.
      E.g. --maskoptions "-species human" has to be written " -species human"
      because "-s..." is taken as the argument -s (I think...)
    - tar command creates myinput.tar.gz. When untarred (tar zxvf myinput.tar.gz) it produces:
          meme/myinput/[subdirs]
      Instead it should give:
          myinput/[subdirs]
      
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-i', '--input', # nargs='?', default=sys.stdout,
                    type= str,
                    required= True,
                    help="""Input bed file containing regions where to search
for motifs. It might be the output of a peak caller which can be parsed by
further arguments
                    """)

parser.add_argument('-o', '--outdir', # nargs='?', default=sys.stdout,
                    type= str,
                    required= False,
                    default= None,
                    help="""Output dir where to send results. Default None will
send output to ./meme/inputfilename/   
                    """)

parser.add_argument('-a', '--awk',
                    type= str,
                    default= '$14 < 5',
                    help="""Filter the input bed using this string passed to awk.
Default '$14 < 5' is suitable for the *.macs_combined.bed output to filter out
peaks with FDR < 5%%
                    """)

parser.add_argument('-s', '--sort',
                    type= str,
                    default= '-k 12,12nr -k 5,5nr',
                    help="""Rank peaks according to this expression passed to
bash/sort. Default is '-k 12,12nr -k 5,5nr'.
Column 12 and 5 are -10*log10(pvalue) and tags at peak summit.
                    """)

parser.add_argument('-n', '--npeaks',
                    type= int,
                    default= 100,
                    help="""Use only the first n regions (after filtering and sorting).
This int passed to bash/head -n
                    """)

parser.add_argument('-e', '--extendby',
                    type= int,
                    default= 100,
                    help="""Extend the coordinates in the input bed file by
these many bases. Default is 100. E.g. if the input bed has the coordinates
of the peak summit "chr1 200 201 ..." they will become "chr1 100 301".
                    """)

parser.add_argument('-r', '--refgenome',
                    type= str,
                    default= '~/reference_seqs/Homo_sapiens_NCBI_v36.3/hsa.fa',
                    help="""Fasta file from which to extract sequences. This
file and the filtered input bed are passed to BEDTools/FastaFromBed 
                    """)

parser.add_argument('-M', '--mask',
                    action= 'store_true',
                    help="""With this option the extracted input sequences are
masked using RepeatMasker. By default, don't maks. 
                    """)

parser.add_argument('--maskoptions',
                    type= unicode,
                    default= '',
                    help="""String of options to pass to RepeatMasker.
Defualt is no options (''). Consider in particular the '-species' option
(e.g. " -species human"). It is necessary *not* to start the string with a
hyphen, start with a whitespace (bug to be fixed).
                    """)

parser.add_argument('-S', '--subsample',
                    type= int,
                    default= 600,
                    help="""Option passed to meme/fasta-subsample: Number of
sequences random sample from the input fasta. The sampled sequences will be
passed to the motif enrichment steps. Default 600 (reduced if this number
exceeds the number of available seqs.)
                    """)

parser.add_argument('--memeopt',
                    type= str,
                    default= '-dna -mod anr -nmotifs 5 -minw 15 -maxw 50 -time 7200 -revcomp',
                    help="""Option passed to meme main program as string.
Do not include the option -o/-oc (output dir) as this is automatically created.

NB: meme default is to use
    * protein alphabet (use -dna for dna)
    * -minw, -maxw     min and max width of the motif: 8, 50
    * -mod zoop as default from web based meme
    * -nmotifs 3
                    
Default: "-dna -mod anr -nmotifs 5 -minw 15 -maxw 50 -time 7200 -revcomp"
""")

parser.add_argument('-m', '--memedb',
                    type= str,
                    default= '/home/berald01/applications/meme/db/motif_databases/JASPAR_CORE_2009.meme',
                    help="""Database of motifs to pass to the meme suite.
Default is JASPAR_CORE_2009 
                    """)

parser.add_argument('-p', '--scp',
                    type= str,
                    default= '',
                    help="""Send (the tar & gzipped) output using scp to this
address and directory. E.g. '143.65.172.17:~/Tritume'. Username is the one
currently logged in. Default is not to do any tar & scp.
                    """)

args = parser.parse_args()

"-------------------------------[ Functions ]--------------------------------- "

def count_fasta_seq(fasta_file):
    """
    Count the number of sequences in the input fasta_file.
    Returns an integer.
    """
    try:
        fh= open(fasta_file)
    except:
        return('chipseq_motif.count_fasta_seq: Cannot open file %s' %(fasta_file))
    ncount= 0
    for n in fh:
        if n.startswith('>'):
            ncount += 1
    fh.close()
    return(ncount)

def extend_bed(fasta_file, extendby= 100, fasta_output= None):
    """
    Extends the coordinates of a bed file by 'extandby' bases.
    Negative coordinates reset to zero.
    
    fasta_file= input bed
    extendby= subtract from 2nd col and add to 3rd column these many bases.
    fasta_output= Output name. If None, overwrite the input file
    """
    fastain= open(fasta_file)
    fastaseq= fastain.readlines()
    fastain.close()
    if fasta_output is None:
        fastaout= open(fasta_file, 'w')
    else:
        fastaout= open(fasta_output, 'w')
    for line in fastaseq:
        line= line.strip('\r\n')
        line= line.split('\t')
        line[1]= str(int(line[1]) - extendby)
        if int(line[1]) < 0:
            line[1]= '0'
        line[2]= str(int(line[2]) + extendby)
        line= '\t'.join(line)
        fastaout.write(line + '\n')
    fastaout.close()
    
" ------------------------[ Clean-up and prepare for outputs ]---------------- "

if args.outdir is None:
    outdir= os.path.splitext(os.path.basename(args.input))[0]
    meme_output_dir= os.path.join('meme', outdir)
else:
    meme_output_dir= args.outdir

print('\nRemoving and creating dir: %s\n' %(meme_output_dir))
shutil.rmtree(meme_output_dir, ignore_errors= True)
os.makedirs(meme_output_dir)

peak_bed=   os.path.join(meme_output_dir, os.path.splitext(os.path.basename(args.input))[0]  + '.top_peaks.bed') ## Where the filtered/sorted input will be
peak_fasta= os.path.join(meme_output_dir, os.path.splitext(os.path.basename(peak_bed))[0] +  '.fa')              ## Where the fasta sequences will be

" ------------------------[ Filter and sort input file]----------------------- "

## The awk|sort and head programs could be executed in the same bash line, but I get the error
## sort: write failed: standard output: Broken pipe
## although the output seems fine.
cmd= """awk '%s' %s | sort %s > %s""" %(args.awk, args.input, args.sort, peak_bed + '.tmp')
print(cmd)
p= subprocess.Popen(cmd, shell= True)
p.wait()
cmd= "head -n %s %s > %s; rm %s" %(args.npeaks, peak_bed + '.tmp', peak_bed, peak_bed + '.tmp')
print(cmd)
p= subprocess.Popen(cmd, shell= True)
p.wait()

if args.extendby > 0:
    extend_bed(peak_bed, args.extendby)

" -----------------------------[ Extract fasta seqs ]------------------------- "
cmd= """fastaFromBed -fi %s -bed %s -fo %s""" %(args.refgenome, peak_bed, peak_fasta)
print(cmd)
p= subprocess.Popen(cmd, shell= True)
p.wait()

" ----------------------------[ Repeat masker ]------------------------------- "

cmd= """RepeatMasker %s %s""" %(args.maskoptions, peak_fasta)
if args.mask:
    print(cmd)
    p= subprocess.Popen(cmd, shell= True)
    p.wait()

" -----------------------------[ Meme module ]-------------------------------- "

if args.mask:
    input_fasta= peak_fasta + '.masked'
else:
    input_fasta= peak_fasta

## Maximum no. seq to analyze
tot_fasta= count_fasta_seq(input_fasta)
if args.subsample < tot_fasta:
    no_subsample= args.subsample
else:
    no_subsample= tot_fasta

print('\nSampling %s out of %s sequences' %(no_subsample, tot_fasta))

cmd= """
fasta-center -len 200 < %(input_fasta)s > %(output_dir)s/seqs-centered;
echo No. sequences in %(output_dir)s/seqs-centered: `grep '>' %(output_dir)s/seqs-centered | wc -l`;
fasta-dinucleotide-shuffle -f %(output_dir)s/seqs-centered -t -dinuc > %(output_dir)s/seqs-shuffled;
cat %(output_dir)s/seqs-centered %(output_dir)s/seqs-shuffled > %(output_dir)s/seqs-centered_w_bg;
fasta-subsample %(output_dir)s/seqs-centered %(no_subsample)s -rest %(output_dir)s/seqs-discarded > %(output_dir)s/seqs-sampled;
""" %{'input_fasta': input_fasta, 'output_dir': meme_output_dir, 'no_subsample': no_subsample}
print(cmd)
p= subprocess.Popen(cmd, shell= True)
p.wait()

cmd= """
meme %(output_dir)s/seqs-sampled -oc %(output_dir)s/meme_out %(meme_opt)s; ## default -minw 6 -maxw 30
tomtom -verbosity 1 -oc %(output_dir)s/meme_tomtom_out -min-overlap 5 -dist pearson -evalue -thresh 0.1 -no-ssc %(output_dir)s/meme_out/meme.txt %(meme_db)s;
dreme2 -v 1 -oc %(output_dir)s/dreme_out -p %(output_dir)s/seqs-centered -n %(output_dir)s/seqs-shuffled -png;
tomtom -verbosity 1 -oc %(output_dir)s/dreme_tomtom_out -min-overlap 5 -dist pearson -evalue -thresh 0.1 -no-ssc %(output_dir)s/dreme_out/dreme.txt %(meme_db)s;
ame --verbose 1 --oc %(output_dir)s/ame_out --fix-partition 800 --bgformat 0 %(output_dir)s/seqs-centered_w_bg %(meme_db)s;
""" %{'input_fasta': input_fasta, 'output_dir': meme_output_dir, 'meme_db': args.memedb, 'meme_opt': args.memeopt}
print(cmd)
p= subprocess.Popen(cmd, shell= True)
p.wait()

" ---------------------[ Compress and ship results ]-------------------------- "

if args.scp != '':
    uname= getpass.getuser()
    cmd= """
tar zcvf %(meme_output_dir)s.tar.gz %(meme_output_dir)s > /dev/null
scp %(meme_output_dir)s.tar.gz %(scp_to)s
""" %{'meme_output_dir': meme_output_dir, 'scp_to': uname + '@' + args.scp}
    print(cmd)
    p= subprocess.Popen(cmd, shell= True)
    p.wait()
    
sys.exit()