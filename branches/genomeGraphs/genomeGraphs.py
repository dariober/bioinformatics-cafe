#!/usr/bin/env python

import argparse
import sys
import os
import tempfile
import shutil
import glob
import multiprocessing
from pycoverage import *
from validate_args import *
import pympileup
import pipes

# MEMO: Add recyclable graphical parameters.
# E.g. we want to add a new col parameter that has to be applied to each track.
# Recycled if not enough items are passed to command line.
# 
# 1. Add arguments to argparse parser. Use nargs= '+' to make it a list:
#   p.add_argument('--col_line', default= ['blue'], nargs= '+', ...)
# 
# 2. Add this arguments to RPlot function. Enclose the list in quoteStringList.
#   This will send the arg to the R script:
# RPlot(..., col_line= quoteStringList(args.col_line), ...)
#
# 3. In R_template.R: Get this argument by assigning to a var and recycle as
#    necessary. col_line is now a vector c('blue', 'red', ...):
# col_line<- recycle(length(inputlist), c(%(col_line)s))

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Produce coverage plots for one or more bam/bedgraph/gtf/bed files at the
    intervals specified in a bed file. 
    
    Plots can be annotated according to a gtf or bed file and decorated with the individual
    nucleotides if a corresponding reference FASTA file is provided.
    
    Intermediate output files, including the R script, can be saved for future inspection.

EXAMPLE:
    ## Plot coverage of all the bam files in current dir in the region(s) in file actb.bed
    ## Annotate plot given a GTF file.
    genomeGraphs.py --ibam *.bam --bed actb.bed

    ## Keep intermediate files:
    genomeGraphs.py --ibam ds05*.bam --bed actb.bed --tmpdir actb

SEE ALSO:
    Documentation at
    http://code.google.com/p/bioinformatics-misc/wiki/coverage_screenshots_docs

    """) # , formatter_class= argparse.RawTextHelpFormatter

# -----------------------------------------------------------------------------
input_args= parser.add_argument_group('Input options', '')

input_args.add_argument('--ibam', '-i',
                   required= False,
                   default= [],
                   nargs= '+',
                   help='''List of input files for which coverage or annotation
should be plotted. Bam files must be sorted and indexed. bedGraph files are
taken as coverage all the other formats (bed, gtf, generic txt) are "annotation".
Input can be gzipped. Metacharacters are expanded by python (`glob`). E.g. to
match all the bam files use '*.bam'. Use '-' to read the list of files from stdin.
                   ''')

input_args.add_argument('--bed', '-b',
                   required= True,
                   help='''Bed or Gtf file with regions to plot optionally gzipped.
Use - to read from stdin.
                   ''')

input_args.add_argument('--slop', '-s',
                   default= ['0.05', '0.05'],
                   nargs= '+',
                   help='''Extend each interval in --bed input. If integer(s),
extend by that many base left and/or right. If float(s), extend by that percent
of interval size (e.g. 0.1 to extend by 10%%). If one value is given, it will be
applied to left and right. If two values, first will be applied to left and second
to right. Must be >= 0.
                   ''')


input_args.add_argument('--fasta', '-f',
                   help='''Fasta file of the reference genome. If given, high
resolution plots will show the reference bases at each position.
                   ''')

input_args.add_argument('--samtools',
                    default= '',
                    help='''Path to samtools. Default is '' which assumes it is
on PATH''')

input_args.add_argument('--parfile', '-pf',
                    default= None,
                    help='''_In prep_: Parameter file to get arguments from.''')

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
                   help='''If given, concatenate all the pdf figures into this
output file (requires PyPDF2 package). 
                   ''')

output_args.add_argument('--tmpdir', '-t',
                    default= None,
                    help='''Directory where to dump temporary files. If not assigned
the tmpdir will be deleted. If set it will *not* be deleted.
                   ''')

output_args.add_argument('--nwinds', '-w',
                   type= int,
                   default= 1000,
                   help='''Maximum number of data-points to plot. If the bed interval
contains more than --nwinds data points it will be divided into equally sized
windows and the counts or bedgraph scores grouped by window. Small value give a coarse
resolution while larger values more jagged profile. Default 1000.
If nwinds < maxseq,  nwinds is reset to maxseq.  
''')

output_args.add_argument('--replot',
                   action= 'store_true',
                   help='''Re-draw plots using the output files from a previous
execution. This option allows to reformat the plots without going thorugh the time
consuming steps of generating pileups, intersections etc.
''')

output_args.add_argument('--rpm',
                   action= 'store_true',
                   help='''Normalize counts by reads per million using library
sizes. Default is to use raw counts. (Only relevant to bam files).
''')

output_args.add_argument('--verbose', '-v',
                   action= 'store_true',
                   help='''Print verbose output. Currently this option only adds
the stdout and stderr from R. Only useful for debugging.
''')

output_args.add_argument('--group_fun',
                   type= str,
                   default= 'mean',
                   help='''The function to group-by windows if the bedgraph or bam files
generate ranges larger than --nwinds. Default 'mean'. See bedtools groupby for
valid options.
''')

# -----------------------------------------------------------------------------
plot_coverage= parser.add_argument_group('Plot of coverage', '')

plot_coverage.add_argument('--maxseq', '-m',
                    default= 100,
                    type= int,
                    help='''The maximum width of the region (bp) to print bases
and to plot each nucleotide in different colour. Default 100 (i.e. regions smaller
than 100 bp will be printed with colour coded nucleotides and the sequence will
be shown).''')

plot_coverage.add_argument('--col_nuc', nargs= '+', default= [''], help='''Colour for the four
nucleotides. _in prep_ File to the colour
code for counts of the ACTGN (counts on + strand) and actgn (counts on - strand).
''')

plot_coverage.add_argument('--no_col_bases', action= 'store_true', help='''Do not colour
individual bases even if the span of the region is < maxseq. I.e. use col_track
and col_track_rev
''')

plot_coverage.add_argument('--col_all', '-c',
                    action= 'store_true',
                    help='''Paint each bar with the colours given in --col_nuc
even if the base matches the reference. Default is to paint only mismatching bases.
Irrelvant if --maxseq is exceeded or a reference is not given to --fasta.
                   ''')

# -----------------------------------------------------------------------------
annotation_args= parser.add_argument_group('Track options',
    'Graphical options for the individual tracks.')

annotation_args.add_argument('--col_line', default= ['blue'], nargs= '+',
    help='''Colour of the profile line of each coverage track. Default blue''')
annotation_args.add_argument('--lwd', default= [0.2], type= float, nargs= '+',
    help='''Line width of the profile line of each coverage track''')
annotation_args.add_argument('--cex', default= 1, type= float,
                             help='''Character expansion for all text. All the other
cex parameters will be based on this''')
annotation_args.add_argument('--col_text_ann', default= 'black', help='''Colour for annotation text (gene names)''')
annotation_args.add_argument('--col_track', default= [''], nargs= '+', help='''Colour for coverage and annotation tracks and for the N base.
Default will assign grey to coverage and firebrick4 to annotation''')
annotation_args.add_argument('--col_track_rev', default= [''], nargs= '+', help='''Relevant to bam files only: Colour for reads on reverse strand. Use NA
to use the same colour as for forward reads (i.e. turn it off)''')

annotation_args.add_argument('--rcode', default= [''], nargs= '+',
    help='''A list of strings, one for each plot and recycled, that will be evaluated as R code at the end of each plot.
Useful to annotate plots. E.g. "abline(h= 10)"''')

# -----------------------------------------------------------------------------
xaxis_args= parser.add_argument_group('Annotation of x-axis',
    'Affect x-axis labelling, range and sequence of interval')

xaxis_args.add_argument('--cex_axis', default= 1, type= float, help='''Character exapansion for the axis annotation.''')
xaxis_args.add_argument('--cex_range', default= 1, type= float, help='''Character exapansion for the range of the x-axis''')
xaxis_args.add_argument('--cex_seq', default= 1, type= float, help='''Character exapansion for the nucleotide sequence''')
xaxis_args.add_argument('--col_seq', default= 'black', help='''Colour for the nucleotide sequence.''')

# -----------------------------------------------------------------------------
plot_layout= parser.add_argument_group('Plot layout', '')

plot_layout.add_argument('--ymax', '-Y',
                    default= ['indiv'],
                    type= str,
                    nargs= '+',
                    help='''Maximum limit of y-axis. Options are:
'indiv': Scale each plot individually to its maximum (default).
'max' all y-axes set to the maximum value of all the coverage plots.
<float>: Set all the plots to this maximum.
'max' cannot be combined with other choices. Float and indiv will be recycled.
                   ''')

plot_layout.add_argument('--ymin', '-y',
                    default= ['min'],
                    type= str,
                    nargs= '+',
                    help='''Minimum limit of y-axis. Options are:
'min' (default) all y-axes set to 0 or the minimum value of all the coverage plots.
'indiv': Scale each plot individually to its minimum.
<float>: Set all the plots to this minimum.''')

plot_layout.add_argument('--ylab', '-yl',
                    default= [''],
                    type= str,
                    nargs= '+',
                    help='''Labels for Y-axis. Recycled.''')

plot_layout.add_argument('--vheights',
                    default= [''],
                    type= str,
                    nargs= '+',
                    help='''List of proportional heights to be passed to R layout(). Recycled.
E.g. if `4 2 1` will make the 1st track twice the size of the 2nd and 4 times the height of 3rd.''')

plot_layout.add_argument('--names', default= None, nargs= '+', help='''List of names for the samples. Default ('') is to use the names of the
bam files with path and .bam extension stripped. Recycled as necessary.''')
plot_layout.add_argument('--cex_names', default= 1, type= float, help='''Character exapansion for the names of the samples''')
plot_layout.add_argument('--col_names', default= ['#0000FF50'], nargs= '+',
    help='''List of colours for the name of each samples. Colours recycled as necessary.
Useful to colour-code samples according to experimemtal design.''')

plot_layout.add_argument('--bg', nargs= '+', default= ['grey95'], help='''List of colours for the plot backgrounds. Recycled as necessary.
Useful to colour code samples sharing the same conditions''')

# -----------------------------------------------------------------------------
figure_size_args= parser.add_argument_group('Global graphical optons',
    'These options affect all the tracks or the figure as a whole')

figure_size_args.add_argument('--title', default= '',
    help= '''Title for each region. Default are the region coords <chrom>:<start>-<end>
You can add a custom title while keeping the coordinates using the syntax --title ':region: My title'
TODO: Make it recyclable/one for each region?''')

figure_size_args.add_argument('--mar', default= [0, 4, 0.2, 1], nargs= 4, type= float, help='''List of 4 floats giving the margins of each plot.''')

figure_size_args.add_argument('--col_grid', default= ['darkgrey'], nargs= '+', help='''Colour for grid(s), recycled.''')

figure_size_args.add_argument('--pwidth', '-W',
                    default= 15,
                    type= float,
                    help='''Width of the plots in cm. Default 15
                    ''')

figure_size_args.add_argument('--pheight', '-H',
                    default= -1,
                    type= float,
                    help='''Height of the figure in cm
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
    if args.parfile:
        pf= read_parfile(args.parfile)
        if not pf:
            sys.exit('Error parsing parameter file %s' %(args.parfile))
        args= assign_parfile(pf, args)
    # Checking arguments
    # ------------------
    if args.ibam == '-' and args.bed == '-':
        sys.exit('stdin passed to *both* --ibam and --bed!')
    try:
        assert validate_ymax(args.ymax)
        assert validate_ymin(args.ymin)
    except AssertionError:
        print('''Invalid arguments passed to ymax (%s) or ymin (%s).''' %(args.ymax, args.ymin))
        sys.exit(1)
    if len(args.slop) > 2:
        print('Only one or two ints or floats must be passed to --slop! Got %s' %(args.slop))
        sys.exit(1)
    if len(args.slop) == 1:
        args.slop.append(args.slop[0])
    slop= []
    try:
        slop.append(int(args.slop[0]))
    except ValueError:
        try:
            slop.append(float(args.slop[0]))
        except:
            sys.exit('Invalid argument passed to --slop. Must be ints or floats: %s' %(args.slop))
    try:
        slop.append(int(args.slop[1]))
    except ValueError:
        try:
            slop.append(float(args.slop[1]))
        except:
            sys.exit('Invalid argument passed to --slop. Must be ints or floats: %s' %(args.slop))
            
    if args.nwinds < args.maxseq:
        print('\nWarning: nwinds (%s) reset to maxseq (%s)' %(args.nwinds, args.maxseq))
        nwinds= args.maxseq
    else:
        nwinds= args.nwinds
    if args.replot and args.tmpdir is None:
        sys.exit('\nCannot replot without a working (--tmpdir) directory!\n')
    if args.ibam == ['-']:
        inputlist_all= [x.strip() for x in sys.stdin.readlines()]
    else:
        inputlist_all= getFileList(args.ibam)
    inputlist= dedupFileList(inputlist_all)
    bamlist= [x for x in inputlist if x.endswith('.bam')]
    nonbamlist= [x for x in inputlist if not x.endswith('.bam')]
    names= parse_names(args.names, inputlist_all)
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
        outdir= os.path.split(args.onefile)[0]
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
        libsizes_dict= pympileup.getLibrarySizes(bamlist, args.samtools)
        libsizes= [libsizes_dict[x] for x in bamlist]
        print(', '.join([str(x) for x in libsizes]))
    if args.bed == '-':
        inbed= sys.stdin
    elif args.bed.endswith('.gz'):
        inbed= gzip.open(args.bed)
    else:
        inbed= open(args.bed)
    inbed= pybedtools.BedTool(inbed).sort() ## inbed is args.bed file handle
    
    # ---------------------[ Pre-filter non-bam files ]-------------------------
    xinbed= pybedtools.BedTool(inbed).each(slopbed, slop)
    nonbam_dict= {}
    for nonbam in nonbamlist:
        nonbam_dict[nonbam]= prefilter_nonbam_multiproc(nonbam= nonbam, inbed= xinbed, tmpdir= tmpdir)
#        proc_list.append({'nonbam':nonbam, 'inbed':inbed, 'tmpdir':tmpdir})
#    pool = multiprocessing.Pool(processes= args.nproc)
#    ori_new= pool.map(prefilter_nonbam_multiproc, proc_list)
#    nonbam_dict= {}
#    for t in ori_new:
#        nonbam_dict[t[0]]= t[1]
    # -----------------------[ Loop thorugh regions ]----------------------------
    for region in inbed:
        print('Processing: %s' %(str(region).strip()))
        bstart= region.start
        bend= region.end
        regname= '_'.join([str(x) for x in [region.chrom, bstart, bend]])
        if region.name != '' and region.name != '.':
            regname = regname + '_' + region.name
        regname= re.sub('[^a-zA-Z0-9_\.\-\+]', '_', regname) ## Get rid of metachar to make sensible file names
        xregion= slopbed(region, slop)
        
        ## --------------------[ Prepare output file names ]-------------------
        
        fasta_seq_name= os.path.join(tmpdir, regname + '.seq.txt')
        if bamlist != []:
            mpileup_name= os.path.join(tmpdir, regname + '.mpileup.bed.txt')
            mpileup_grp_name= os.path.join(tmpdir, regname + '.grp.bed.txt')
        else:
            mpileup_name= ''
            mpileup_grp_name= ''
        if nonbamlist != []:
            non_bam_name= os.path.join(tmpdir, regname + '.nonbam.bed.txt')
        else:
            non_bam_name= ''
        pdffile= os.path.join(tmpdir, regname + '.pdf')
        final_pdffile= os.path.join(outdir, regname + '.pdf')
        rscript= os.path.join(tmpdir, regname + '.R')
        if not args.replot:
            prepare_reference_fasta(fasta_seq_name, args.maxseq, xregion, args.fasta) ## Create reference file even if header only
            ## ----------------------- BAM FILES -------------------------------
            ## At the end of this session you have *.grp.bed.txt (matrix-like
            ## file read by R)
            regionWindows= False ## bed interval divided into nwinds intervals by bedtools windowMaker
            if bamlist != []:
                if not regionWindows:
                    regionWindows= makeWindows(xregion, nwinds)
                pympileup.bamlist_to_mpileup(mpileup_name, mpileup_grp_name, bamlist, xregion, nwinds, args.fasta, args.rpm, regionWindows, samtools= args.samtools, groupFun= args.group_fun) ## Produce mpileup matrix
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
                non_bam_fh= open(non_bam_name, 'w') ## Here all the files concatenated.
                for x in nonbamlist:
                    nonbam= nonbam_dict[x]
                    if x.endswith('.bedGraph') or x.endswith('.bedGraph.gz'):
                        """Bedgraph needs to go to tmp file because you don't know if
                        it has to be compressed by windows or not. <- This should be
                        changed: You have already intersected the nonbam files with
                        the bed regions.
                        """
                        tmpfh= tempfile.NamedTemporaryFile(dir= tmpdir, suffix= 'nonbam.tmp.bed', delete= False)
                        tmp_name= tmpfh.name
                        #tmpfh= open(tmp_name, 'w')
                        nlines= prepare_nonbam_file(nonbam, tmpfh, xregion, use_file_name= x) ## Write to fh the overlaps btw nonbam and region. Return no. lines
                        tmpfh.close()
                        if nlines > nwinds:
                            if not regionWindows:
                                regionWindows= makeWindows(xregion, nwinds)
                            compressBedGraph(regionWindows, tmp_name, use_file_name= x, bedgraph_grp_fh= non_bam_fh, col_idx= 4 + len(pympileup.COUNT_HEADER),
                                groupFun= args.group_fun)
                        else:
                            fh= open(tmp_name)
                            for line in fh:
                                non_bam_fh.write(line)
                        os.remove(tmp_name)
                    else:
                        nlines= prepare_nonbam_file(nonbam, non_bam_fh, xregion, use_file_name= x) ## Write to fh the overlaps btw nonbam and region. Return no. lines
                non_bam_fh.close()
            else:
                non_bam_name= ''
        # ----------------------------------------------------------------------
        # Plotting 
        # ----------------------------------------------------------------------
        outputPDF.append(pdffile)
        rgraph= RPlot(
              inputlist= quoteStringList(inputlist_all),
              count_header= quoteStringList(pympileup.COUNT_HEADER),
              pdffile= pdffile,
              rscript= rscript,
              mcov= mpileup_grp_name,
              nonbam= non_bam_name,
              refbases= fasta_seq_name,
              title= args.title,
              pheight= args.pheight,
              pwidth= args.pwidth,
              psize= args.psize,
              ylab= quoteStringList(args.ylab),
              bstart= bstart,
              bend= bend,
              xlim1= xregion.start,
              xlim2= xregion.end,
              maxseq= args.maxseq,
              ymax= quoteStringList(args.ymax),
              ymin= quoteStringList(args.ymin),
              chrom= xregion.chrom,
              vheights= quoteStringList(args.vheights),
              cex= args.cex,
              cex_axis= args.cex_axis,
              col_line= quoteStringList(args.col_line),
              lwd= quoteStringList(args.lwd),
              col_track= quoteStringList(args.col_track),
              col_track_rev= quoteStringList(args.col_track_rev),
              col_nuc= quoteStringList(args.col_nuc),
              no_col_bases= args.no_col_bases,
              bg= quoteStringList(args.bg),
              col_grid= quoteStringList(args.col_grid),
              col_text_ann= args.col_text_ann,
              names= quoteStringList(names),
              col_names= quoteStringList(args.col_names),
              cex_names= args.cex_names,
              cex_range= args.cex_range,
              cex_seq= args.cex_seq,
              col_seq= args.col_seq,
              mar= ', '.join([str(x) for x in args.mar]),
              col_all= args.col_all,
              rcode= quoteStringList(args.rcode)
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
            shutil.copyfile(pdffile, final_pdffile)
    if onefile:
        catPdf(in_pdf= outputPDF, out_pdf= args.onefile)
    for f in nonbam_dict:
        os.remove(nonbam_dict[f])
    if args.tmpdir is None:
        shutil.rmtree(tmpdir)
if __name__ == '__main__':
    main()
    sys.exit()
