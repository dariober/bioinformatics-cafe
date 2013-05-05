#!/usr/bin/env python

import argparse
import sys
import os
import tempfile
import subprocess
import shutil
from pycoverage import *

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Produce coverage plots for one or more bam files at the positons specified in
    a bed file. Plots are written in pdf format, one per region or concatenated in
    a single file.
    
    Plots can be annotated according to a GTF file and decorated with the individual
    nucleotides if a corresponding refernce FASTA file is provided.
    
    Intermediate output files, including the R script, can be saved for future inspection.

EXAMPLE:
    ## Plot coverage of all the bam files in current dir in the region(s) in file actb.bed
    ## Annotate plot given a GTF file.
    coverage_screenshots.py --ibam *.bam --gtf genes.gtf.gz --bed actb.bed

    ## Keep intermediate files:
    coverage_screenshots.py --ibam ds05*.bam --gtf genes.gtf.gz --bed actb.bed --tmpdir actb

SEE ALSO:
    Documentation at
    http://code.google.com/p/bioinformatics-misc/wiki/coverage_screenshots_docs

    """, formatter_class= argparse.RawTextHelpFormatter)

# -----------------------------------------------------------------------------
input_args= parser.add_argument_group('Input options', '')

input_args.add_argument('--ibam', '-i',
                   required= False,
                   default= [],
                   nargs= '+',
                   help='''List of bam files, sorted and indexed to visualize.
Metacharacters are expanded by python (`glob`). E.g. to match all the bam files
use '*.bam'. 
                   ''')

input_args.add_argument('--bed', '-b',
                   required= True,
                   help='''Bed file with regions to plot. 
                   ''')

input_args.add_argument('--fasta', '-f',
                   help='''Fasta file of the reference genome. If given, high
resolution plots will show the reference bases at each position.
                   ''')

input_args.add_argument('--gtf', '-g',
                    help='''GTF file to fetch annotation from. E.g. gene.gtf in
iGenomes. Can be gzip compressed.
                   ''')

input_args.add_argument('--samtools',
                    default= '',
                    help='''Path to samtools. Default is '' which assumes it is
on PATH''')

# -----------------------------------------------------------------------------
output_args = parser.add_argument_group('Output options', '')

output_args.add_argument('--outdir', '-d',
                   required= False,
                   default= None,
                   help='''Output directory for the pdf files. Default to current
dir. NB: Not to be confused with --tmpdir where temp file go.
                   ''')

output_args.add_argument('--onefile', '-o',
                   required= False,
                   default= None,
                   help='''Concatenate the output PDF files into a single one
passed as argument (requires PyPDF2 package). 
                   ''')

output_args.add_argument('--tmpdir', '-t',
                    default= None,
                    help='''Directory where to dump temporary files. By default
python will find one which will be deleted at the end of the execution. If set
it will *not* be deleted.
                   ''')

output_args.add_argument('--replot',
                   action= 'store_true',
                   help='''Re-use the output files from a previous execution. Just
redraw the plots using different graphical parameters.
This option allows to reformat the plots without going thorugh the time consuming
steps of generating pileups etc.
''')

output_args.add_argument('--rpm',
                   action= 'store_true',
                   help='''Normalize counts by reads per million using library
sizes. Default is to use raw counts. 
''')

output_args.add_argument('--verbose', '-v',
                   action= 'store_true',
                   help='''Print verbose output. Currently this option only adds
the stdout and stderr from R. Only useful for debugging.
''')

output_args.add_argument('--nwinds', '-w',
                   type= int,
                   default= 1000,
                   help='''Maximum number of data-points to plot. If the bed interval
is larger than --nwinds it will be divided into equally sized windows and counts
averaged by window. Small value give a coarse resolution while larger values more
jagged profile. Default 1000. If nwinds < maxres,  nwinds is reset to maxres.  
''')

# -----------------------------------------------------------------------------
plot_coverage= parser.add_argument_group('Plot of coverage', '')

plot_coverage.add_argument('--col_nuc', nargs= 4, default= '', help='''List of 4 R colours for bars of A, C, G, and T.''')
plot_coverage.add_argument('--col_cov', default= 'grey', help='''Colour for coverage of pooled bases and for N''')
plot_coverage.add_argument('--col_all', '-c',
                    action= 'store_true',
                    help='''Paint each bar with the colours given in --col_nuc
even if the base matches the reference. Default is to paint only mismatching bases.
Irrelvant if --maxres is exceeded or a reference genome is not provided with --fasta.
                   ''')

plot_coverage.add_argument('--maxres', '-m',
                    default= 100,
                    type= int,
                    help='''The maximum width of the region (bp) to print bases
and to plot each nucleotide in different colour. Default 100 (i.e. regions smaller
than 100 bp will be printed with colour coded nucleotides and the sequence will
be shown).''')

# -----------------------------------------------------------------------------
annotation_args= parser.add_argument_group('Format of annotation', '')

annotation_args.add_argument('--col_text_ann', default= 'black', help='''Colour for annotation text (gene names)''')
annotation_args.add_argument('--col_ann', default= 'firebrick4', help='''Colour for annotation bars (exons, CDS etc.)''')
annotation_args.add_argument('--names', default= '', nargs= '+', help='''List of names for the samples. Default ('') is to use the names of the
bam files with path and .bam extension stripped. Recycled as necessary.''')
annotation_args.add_argument('--cex_names', default= 0.9, type= float, help='''Character exapansion for the names of the samples''')
annotation_args.add_argument('--col_names', default= ['#0000FF50'], nargs= '+',
    help='''List of colours for the name of each samples. Colours recycled as necessary.
Useful to colour-code samples according to experimemtal design.''')
annotation_args.add_argument('--cex_range', default= -1, type= float, help='''Character exapansion for the text range of plot''')
annotation_args.add_argument('--line_range', default= 2, type= float, help='''Distance of range bar from x-axis. In R's line units''')
annotation_args.add_argument('--cex_seq', default= -1, type= float, help='''Character exapansion for the nucleotide sequence''')
annotation_args.add_argument('--line_seq', default= 3.5, type= float, help='''Distance of nucleotide sequence bar from x-axis. In R's line units''')
annotation_args.add_argument('--col_seq', default= 'black', help='''Colour for the nucleotide sequence.''')

# -----------------------------------------------------------------------------
plot_layout= parser.add_argument_group('Plot layout', '')

plot_layout.add_argument('--ylim', '-y',
                    default= 'max',
                    type= str,
                    help='''How the maximum value for the y-axis should be set.
The lower limit of the y-axis is always 0. Options are:
'max' (default) all y-axes set to the maximum value of all the plots.
'indiv': Scale each plot individually to its maximum.
<float>: Set all the plots to this maximum (E.g. 1000).
                   ''')

plot_layout.add_argument('--cex_axis', default= -1, type= float, help='''Character exapansion for the axis annotation (cex.axis in R).
Use negative value to set default.''')

plot_layout.add_argument('--bg', nargs= '+', default= ['grey95'], help='''List of colours for the plot backgrounds. Recycled as necessary.
Useful to colour code samples sharing the same conditions''')
plot_layout.add_argument('--nogrid', action= 'store_true', help='''Do not plot grid''')
plot_layout.add_argument('--oma', default= [5, 1.1, 3, 1.1], nargs= 4, type= float, help='''List of 4 floats giving the outer margins of the plot.
Default 4 1.1 3 1.1''')
plot_layout.add_argument('--mar', default= [0.5, 4, 0.5, 1], nargs= 4, type= float, help='''List of 4 floats giving the margins of each plot.
Default 0.5 4 0.5 1''')

# -----------------------------------------------------------------------------
figure_size_args= parser.add_argument_group('Figure and font size', '')

figure_size_args.add_argument('--pheight', '-H',
                    default= -1,
                    type= float,
                    help='''Height of *each* plot in cm. Default 6
                   ''')

figure_size_args.add_argument('--pwidth', '-W',
                    default= 15,
                    type= float,
                    help='''Width of the plots in cm. Default 24
                   ''')

figure_size_args.add_argument('--psize', '-p',
                    default= 10,
                    type= float,
                    help='''Pointsize for R pdf() function. Sizes
between 9 and 12 should suite most cases. Default 10.
                   ''')

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def main():
    args = parser.parse_args()

    if args.nwinds < args.maxres:
        nwinds= args.maxres
    else:
        nwinds= args.nwinds
    if args.replot and args.tmpdir is None:
        sys.exit('\nCannot replot without a working (--tmpdir) directory!\n')
    inputlist= getFileList(args.ibam)
    bamlist= [x for x in inputlist if x.endswith('.bam')]
    nonbamlist= [x for x in inputlist if not x.endswith('.bam')]
    if not args.replot:
        print('\nFiles to analyze (%s found):\n%s\n' %(len(inputlist), ', '.join(inputlist)))
    if len(inputlist) == 0 and not args.replot:
        sys.exit('No file found!\n')
    # Output settings
    # ---------------
    if (args.outdir is not None) and (args.onefile is not None):
        sys.exit('''\nSpecify either --outdir (for one file for each bed region) OR
--onefile (for one single concatenated file).\n''')
    
    onefile= False
    if args.outdir is None and (args.onefile is None):
        outdir= os.getcwd()
    elif args.onefile is None:
        outdir= args.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    else:
        onefile= True
        try:
            import PyPDF2
        except ImportError:
            sys.exit('''\nModule PyPDF2 could not be imported. Eiher installed it
(see https://pypi.python.org/pypi/PyPDF2 ) or avoid using the --onefile option.\n''')
            
    ## Temp dir to dump intermediate files.
    ## -------------------------------------------------------------------------
    if args.tmpdir is None:
        tmpdir= tempfile.mkdtemp(suffix= '_coverageViewer')
    else:
        tmpdir= args.tmpdir
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
    
    outputPDF= [] ## List of all the pdf files generated. Used only for --onefile
        
    ## -------------------------------------------------------------------------
    if args.rpm and not args.replot:
        sys.stdout.write('Getting library sizes... ')
        libsizes_dict= getLibrarySizes(bamlist, args.samtools)
        libsizes= [libsizes_dict[x] for x in bamlist]
        print(', '.join([str(x) for x in libsizes]))
    inbed= pybedtools.BedTool(args.bed)
    for region in inbed:
        print('Processing: %s' %(str(region).strip()))
        regname= '_'.join([str(x) for x in [region.chrom, region.start, region.end]])
        if region.name != '':
            regname = regname + '_' + region.name 
        ## --------------------[ Prepare output file names ]-------------------
        fasta_seq_name= os.path.join(tmpdir, regname + '.seq.txt')
        if args.gtf:
            annot_file= os.path.join(tmpdir, regname + '.annot.txt')
        else:
            annot_file= ''
        mpileup_name= os.path.join(tmpdir, regname) + '.mpileup.bed.txt'
        mpileup_grp_name= os.path.join(tmpdir, regname) + '.grp.bed.txt'
        non_bam_name= os.path.join(tmpdir, regname) + '.nonbam.bed.txt'
        pdffile= os.path.join(tmpdir, regname + '.pdf')
        rscript= os.path.join(tmpdir, regname + '.R')
        if not args.replot:
            prepare_reference_fasta(fasta_seq_name, args.maxres, region, args.fasta) ## Create reference file even if header only
            ## ----------------------- BAM FILES -------------------------------
            ## At the end of this session you have *.grp.bed.txt (matrix-like
            ## file read by R)
            if bamlist != []:
                bamlist_to_mpileup(mpileup_name, mpileup_grp_name, bamlist, region, args.fasta, args.rpm, nwinds, samtools= args.samtools) ## Produce mpileup matrix
            else:
                mpileup_grp_name= ''
            ## ----------------------NON BAM FILES -----------------------------
            ## Produce coverage and annotation files for non-bam files. One output
            ## file prooduced with format
            ## chrom, start, end, file_name, A, C, G, T, Z.
            ## NB: A,C,G,T are always NA. We keep them only for compatibility
            ## with the output form BAM files. The `score` or `name` column from
            ## bed files (4th) goes to column Z.
            ## file_name has the name of the file as it has been passed to --ibam
            if nonbamlist != []:
                non_bam_fh= open(non_bam_name, 'w')
                non_bam_fh.write('\t'.join(['chrom', 'start', 'end', 'file_name', 'A', 'C', 'G', 'T', 'Z', 'feature', 'name', 'strand']) + '\n')
                for nonbam in nonbamlist:
                    prepare_nonbam_file(nonbam, non_bam_fh, region)
                non_bam_fh.close()
            else:
                non_bam_name= ''
        # ----------------------------------------------------------------------
        # Plotting 
        # ----------------------------------------------------------------------
        outputPDF.append(pdffile)
        if args.rpm:
            ylab= 'Reads per million'
        else:
            ylab= 'Read count'
        ## Memo: All the args passed to RPlot() become part of the R script.
        rgraph= RPlot(
              inputlist= quoteStringList(args.ibam),
              pdffile= pdffile,
              rscript= rscript,
              plotname= regname,
              mcov= mpileup_grp_name,
              nonbam= non_bam_name,
              refbases= fasta_seq_name,
              pheight= args.pheight,
              pwidth= args.pwidth,
              psize= args.psize,
              ylab= ylab,
              xlim1= region.start,
              xlim2= region.end,
              gtf= annot_file,
              maxres= args.maxres,
              ylim= args.ylim,
              cex_axis= args.cex_axis,
              col_cov= args.col_cov,
              col_nuc= quoteStringList(args.col_nuc),
              bg= quoteStringList(args.bg),
              nogrid= args.nogrid,
              col_text_ann= args.col_text_ann,
              col_ann= args.col_ann,
              names= quoteStringList(args.names),
              col_names= quoteStringList(args.col_names),
              cex_names= args.cex_names,
              cex_range= args.cex_range,
              cex_seq= args.cex_seq,
              line_range= args.line_range,
              line_seq= args.line_seq,
              col_seq= args.col_seq,
              oma= ', '.join([str(x) for x in args.oma]),
              mar= ', '.join([str(x) for x in args.mar]),
              col_all= args.col_all
              )
        
        
        if rgraph['stderr'] != '':
            print('\nExpection in executing R script "%s"\n' %(rscript))
            print(rgraph['stdout'])
            sys.exit(rgraph['stderr'])
        if args.verbose:
            print(rgraph['stderr'])
            print(rgraph['stdout'])
        if not onefile and tmpdir != outdir:
            ## Copy PDFs from temp dir to output dir. Unless you want them in onefile or
            ## if the final destination dir has been set to be also the tempdir
            shutil.copyfile(pdffile, os.path.join(outdir, regname + '.pdf'))
    if onefile:
        catPdf(in_pdf= outputPDF, out_pdf= args.onefile)
    if args.tmpdir is None:
        shutil.rmtree(tmpdir)
if __name__ == '__main__':
    main()
    sys.exit()
