import sys
import pybedtools
import os
import tempfile
import subprocess
import glob
import gzip
import re
import csv
import pympileup
import copy
import genomeGraphs
import pycoverage

def read_parfile(parfile):
    """Read the paramater file and return a dictionary with {'param': args}.
    """
    allowed_args= sorted(['ibam',
                   'col_line',
                   'lwd',
                   'col_track',
                   'col_track_rev',
                   'col_text_ann',
                   'ymax',
                   'ymin',
                   'ylab',
                   'cex_lab',
                   'vheights',
                   'names',
                   'cex_names',
                   'col_names',
                   'col_grid',
                   'col_mark',
                   'bg',
                   'rcode',
                   'overplot'])
    fin= open(parfile, "rU")
    fargs= csv.DictReader(fin, delimiter= '\t')
    fDict= {}

    for k in fargs.fieldnames:
        if k not in allowed_args:
            print('\nParameter "%s" cannot be passed by the parameter file.\nAllowed parameters are:\n    %s\n' %(k, '\n    '.join(allowed_args)))
            return(False)
    for line in fargs:
        for k in fargs.fieldnames:
            v= line[k]
            if k in fDict:
                fDict[k].append(v)
            else:
                fDict[k]= [v]
    ## Trim lists:
    for k in fDict:
        v= fDict[k]
        for i in range(len(v))[::-1]:
            if v[i] is None or v[i] == '':
                del v[i]
            else:
                break
    return(fDict)

def stdin_inbed_to_fh(stdin):
    """Parse the lines in fh stdin to write them in bed format to tempfile.
    Returns file tempfile.
    Allowed input in stdin:
    Tab separated "chr\tstart\tend\tother-fields"
    Space sperated "chr start end other-fields"
    """
    bedlines=  stdin.readlines()
    bedlines= [x.strip().split() for x in bedlines if x.strip() != '']
    tmpf= tempfile.NamedTemporaryFile(delete= False)
    for line in bedlines:
        tmpf.write('\t'.join(line[0:4]) + '\n')
    tmpf.close()
    return(tmpf)

class SlopError(Exception):
    pass

def slopbed(interval, slop):
    """Extend bed interval by given slop.
    interval:
        pybedtool.Interval to extend or list or tuple with first three items:
        <str chrom> <int start> <int end>
    slop:
        A list or tuple of two integers or floats. Integers
        will expand the interval by that many bases while floats will expand as
        percentages of feature size (e.g. 0.1 to expand by 10 percent). First
        value applied to left and second to right. Must be 0 or positive.
    Return:
        Same list or pybedtool interval as in input with coordinates extended.

    interval = iter(pybedtools.BedTool('chr1 1 100 asdf 0 + a b c d', from_string=True)).next()

    """
    if type(interval) == pybedtools.cbedtools.Interval:
        left= interval.start
        right= interval.end
    elif type(interval) in (list, tuple):
        left= interval[1]
        right= interval[2]
    else:
        raise SlopError('Invalid type of interval. Expected pybedtools.cbedtools.Interval or list or tuple. Got %s' %(type(interval)))
    if slop[0] < 0 or slop[1] < 0:
        raise SlopError('Slop must be positive int or float. Got %s' %(slop))
    
    size= right - left

    ## Left feature    
    if type(slop[0]) == int:
        xleft= left - slop[0]
    elif type(slop[0]) == float:
        xleft= round(left - (size * slop[0]))
    else:
        raise SlopError('Slop must be of type int or float. Got %s' %(type(slop[0])))
    if xleft < 0:
        xleft= 0
    ## Right feature
    if type(slop[1]) == int:
        xright= right + slop[1]
    elif type(slop[1]) == float:
        xright= round(right + (size * slop[1]))
    else:
        raise SlopError('Slop must be of type int or float. Got %s' %(type(slop[1])))

    if type(interval) == pybedtools.cbedtools.Interval:
        xinterval= copy.copy(interval)
        xinterval.start= xleft
        xinterval.end= xright
    else:
        xinterval= [x for x in interval]
        xinterval[1]= xleft
        xinterval[2]= xright
    return(xinterval)

def assign_parfile(pardict, args):
    """Assign to the parser object args the arguments in dictionary pardict.
    Return:
        Parser object args with parameters updated.
    """
    for k in pardict:
        args.__dict__[k]= pardict[k]
    return(args)
        
def getFileList(files):
    """Expand the list of files using glob. Return a list of unique files.
    """
    inputlist_dup= []
    for bam in files:
        inputlist_dup.extend(glob.glob(bam))
    return(inputlist_dup)

def dedupFileList(x):
    """Remove duplicates from list x
    """
    inputlist= []
    for a in x:
        if a not in inputlist:
            inputlist.append(a)
    return(inputlist)

def prefilter_nonbam_multiproc(inbed, nonbam, tmpdir, sorted):
    """For each filename in nonbamlist do the intersection with the regions in
    inbed. Produce filtered files to tmpdir and return a dict of original filenames
    and the filtered name
    inbed:
        pybedtools.BedTool() object of regions to intersect with the non bam files
    nonbamlist:
        List of non-bam files
    tmpdir:
        Where intersected files will go. Output name is
        tmpdir/x_<original name w/o .gz>
    sorted:
        True/False passed to intersectbed(sorted) to use chromsweep algorithm
    Return:
        Name of the (temp) file with the filetered regions.
#        dict of original names:filtered names
#        E.g. {'bedgraph/profile.bedGraph.gz': 'wdir/x_profile.bedGraph', 'annotation/genes.gtf.gz': 'wdir/x_genes.gtf'}
#        Output files are *sorted*
    """
    bname= re.sub('\.gz$', '', os.path.split(nonbam)[1])
    fn= tempfile.NamedTemporaryFile(dir= tmpdir, suffix= '_' + bname, delete= False)
    x_name= fn.name
    pynonbam= pybedtools.BedTool(nonbam)
    nonbam_x_inbed= pybedtools.BedTool().intersect(a= pynonbam, b= inbed, u= True, sorted= sorted)
    if nonbam_x_inbed.count() == 0:
        pass        
    else:
        nonbam_x_inbed= nonbam_x_inbed.sort().saveas(x_name)
    #ori_new= (nonbam, x_name)
    return(x_name)


def makeWindows(region, n):
    """Divide a region in n windows. If the region size is < n, than each postion
    is returned.
    region:
        pybedtools interval object or string as '<chrom>\t<start>\t<end>'
    n:
        Number of windows to divide region into
    fn:
        File name to write to
    Return:
        Output of pybedtools.BedTool().window_maker
    """
    tmp= tempfile.NamedTemporaryFile(delete= False, suffix= '.bed.txt', prefix= 'windowMaker_')
    tmp.write(str(region))
    tmp.close()
    regionWinds= pybedtools.BedTool().window_maker(b= tmp.name, n= n, stream= False)
    os.remove(tmp.name)
    return(regionWinds)

def mergePDF(filenames, output_filename):
    output = PyPDF2.PdfFileWriter()
    for filename in filenames:
        input = PyPDF2.PdfFileReader(file(filename, "rb"))
        for page in input.pages:
            output.addPage(page)
    outputstream = file(output_filename, "wb")
    output.write(outputstream)
    outputstream.close()

def getRefSequence(fasta, region):
    """Read the fasta file and extract the region given in bedtool interval region.
    Return:
        List of tuples with inner tuple ['chrom', 'start', 'end', 'base']
    NB: You need to reduce the start by 1 because fastaFromBed seems to be 1-based.
    """
    region_str= '\t'.join([region.chrom, str(region.start), str(region.end)])
    bedregion= pybedtools.BedTool(str(region_str), from_string= True)
    seq = bedregion.sequence(fi=fasta, tab= True)
    seq= open(seq.seqfn).read().split('\t')
    seqstring= seq[1].strip()
    seq_table= zip(
        [region.chrom] * len(seqstring), ## Column of chrom
        range(region.start, region.end + 1), ## Column of start pos
        range(region.start + 1, region.end + 2), ## Col of end pos
        list(seqstring) ## Col of bases
        )
    return(seq_table)

def quoteStringList(x):
    '''Join the list of string x in a string where each element is double quoted and
    separted by comma & space.
    Useful to convert python's list of strings to R character vectors
    E.g. x= ['blue', 'black']
    quoteStringList(x) >>> '"blue", "black"'
    NB:
        Numbers are converted to strings and quoted as well! None is converted to
        'NA'.
    NB2:
        This function is not meant to cope well with strings containing double quotes,
        weird metachars etc.
    '''
    s= ''
    for y in x:
        if y is None:
            y= 'NA'
        s= s + '"' + str(y) + '", '
    s= s.strip(', ')
    return(s)

def RPlot(**kwargs):
    """Write to file the R script to produce the plots and execute it using Rscript
    kwargs: 
        Arguments that will be interpolated in the string that make up the script.
    Return:
        A dictionary {'stdout': 'stderr':} with the captured output from Rscript.
    """
    rout= open(kwargs['rscript'], 'w')
    rin= os.path.join(os.path.split(pycoverage.__file__)[0], 'R_template.R')
    rtemplate= open(rin).read()
    rplot= rtemplate %kwargs
    rout.write(rplot)
    rout.close()
    p= subprocess.Popen('Rscript %s' %(kwargs['rscript']), stdout= subprocess.PIPE, stderr= subprocess.PIPE, shell= True)
    stdout, stderr= p.communicate()
    if stderr != '':
        print(stderr)
    return({'stdout':stdout, 'stderr': stderr})
        
def catPdf(in_pdf, out_pdf):
    """Concatenate the PDF files in list `in_pdf` into the single file `out_pdf`:
    Return:
        True if successful False otherwise.
    See also:
        http://www.blog.pythonlibrary.org/2010/05/15/manipulating-pdfs-with-python-and-pypdf/
    """
    import PyPDF2
    output = PyPDF2.PdfFileWriter()
    for pdf in in_pdf:
        pdfOne = PyPDF2.PdfFileReader(file(pdf, "rb"))
        output.addPage(pdfOne.getPage(0))
    outputStream = file(out_pdf, "wb")
    output.write(outputStream)
    outputStream.close()

def prepare_nonbam_file(infile_name, outfile_handle, region, use_file_name):
    """Intersect a bed interval with the bed, bed-like or gtf file. Each intersected
    line is written to outfile_handle in the format:
        chrom start end file_name A C G T Z feature strand
    Where A, C, G, T are always NA. Z is the value in column 4 of input, if
    numeric or NA if string (e.g. from gtf annotation). <feature> is NA if Z is numeric
    or CDS or exon otherwise. 
     
    Output an annotation file suitable for R plotting.
    infile_name:
        Annotation or covergae file. Can be gtf (annotation) or any 4-column
        tab delimited file with chrom, start, end, score/name.
    region:
        A pybedtools.Interval to intersect to the gtf. Intersected feature will be
        sent to output
    outfile_handle:
        Open file handle to write to.
    use_file_name:
        Use this string as filename. Must match the original name passed to --ibam 
    Return:
        Number of intervals that overlap `region` (n lines)
    """
    infile= pybedtools.BedTool(infile_name)
    nlines= 0

    if infile.head(n= 1, as_string= True) == '' :
        """If there are no intersecting features return 0
        """
        pass
    else:
        region_x_infile= pybedtools.BedTool().intersect(a= infile, b= pybedtools.BedTool(str(region), from_string= True), sorted= True)
        for line in region_x_infile:
            if line.name == '':
                name= 'NA'
            else:
                name= line.name
            if line.strand == '':
                strand= '.'
            else:
                strand= line.strand
            if infile_name.endswith('.gtf') or infile_name.endswith('.gtf.gz'):
                outline= [line.chrom, line.start - 1, line.end, use_file_name] + ['NA'] * len(pympileup.COUNT_HEADER) + [line.fields[2], name, strand] ## NA for padding depthACTGNZactgnz
            elif infile_name.lower().endswith('.bedgraph') or infile_name.lower().endswith('.bedgraph.gz'):
                outline= [line.chrom, line.start, line.end, use_file_name] + ['0'] * (len(pympileup.COUNT_HEADER)-1) + [line.name, 'coverage', 'NA', strand] ## The column to plot is line.name, ['0']*9 is for padding
            else:
                outline= [line.chrom, line.start, line.end, use_file_name, ] + ['NA'] * len(pympileup.COUNT_HEADER) + ['generic', line.name, strand]
            outfile_handle.write('\t'.join([str(x) for x in outline]) + '\n')
            nlines += 1
    return(nlines)

def prepare_reference_fasta(fasta_seq_name, maxseq, region, fasta):
    """Output a reference file with the base at each postion in
    the region interval
    fasta_seq_name:
        Output name for reference file. Will have format <chrom> <pos> <refbase>
        included header.
    maxseq:
        Maximum size of the interval to extract sequence. If exceeded, only the header
        will be in ouput file.
    region:
        pybedtools region. All the bases in this interval will be sent to
        fasta_seq_name
    fasta:
        Refernce FASTA file from where to extract sequence
    Returns:
        True on success. Side effect produce the reference file *.seq.txt
    """
    region_seq= open(fasta_seq_name, 'w')
    region_seq.write('\t'.join(['chrom', 'start', 'end', 'base']) + '\n')
    if ((region.end - region.start) <= maxseq) and fasta:
        fasta_seq= getRefSequence(fasta, region)
        for line in fasta_seq:
            region_seq.write('\t'.join([str(x) for x in line]) + '\n')
    region_seq.close()
    return(True)
    
def compressBedGraph(regionWindows, bedgraph_name, use_file_name, bedgraph_grp_fh, col_idx= 4 + len(pympileup.COUNT_HEADER), groupFun= 'mean'):
    """Compress a bedgraph by dividing it in n windows and averaging windows
    regionWindows:
        bed interval divided into n windows by pycoverage.makeWindows
    bedgraph_name:
        bedgraph file to compress and name for output compressed bedgraph.
    use_file_name:
        Put this file name in the output line. Must be the same as original
        input.
    bedgraph_grp_fh:
        Output file handle to write to
    col_idx:
        Column index 1 BASED to summarize.
    groupFun:
        Apply this function to group-by windows. This opt passed to bedtools
        groupby. Check there for valid options

    """
    ## Assign to each pileup position its window --------------------------
    if bedgraph_name.endswith('.gz'):
        ncols= len(gzip.open(bedgraph_name).readline().split('\t')) ## get number of columns in this bedgraph
    else:
        ncols= len(open(bedgraph_name).readline().split('\t')) ## get number of columns in this bedgraph
    bedgraph_winds= pybedtools.BedTool(bedgraph_name).intersect(regionWindows, wb= True)
    ## Aggregate counts in each position by window: -----------------------
    wind_idx= [ncols+1, ncols+2, ncols+3] ## These are the indexes of the columns containing the windows. 1 based!
    bedgraph_grp= bedgraph_winds.groupby(g= wind_idx, c= col_idx, o= groupFun, stream= False)
    for line in bedgraph_grp:
        outline= [line.chrom, str(line.start), str(line.end), use_file_name] +  ['0'] * (len(pympileup.COUNT_HEADER)-1) + [line.name, 'coverage', 'NA', '.']
        bedgraph_grp_fh.write('\t'.join(outline) + '\n')

#def format_paramater_dict(pardict, parser):
#    """Format the arguments in dict pardict, tytpically returned by read_parfile(),
#    to comform to the parser.
#    pardict:
#        Dictionary of {parameter:argument}
#    parser:
#        Parser obejct from argparse
#    """
#    sys.exit()
#    argsDict= parser.__dict__['_option_string_actions']
#    for k in argsDict:
#        store_action= argsDict[k].__dict__
#        nargs= store_action
#        print(k, nargs)

#def prefilter_nonbam(inbed, nonbamlist, tmpdir):
#    """For each filename in nonbamlist do the intersection with the regions in
#    inbed. Produce filtered files to tmpdir and return a dict of original filenames
#   and the filtered name
#    inbed:
#        pybedtools.BedTool() object of regions to intersect with the non bam files
#    nonbamlist:
#        List of non-bam files
#    tmpdir:
#        Where intersected files will go. Output name is
#        tmpdir/x_<original name w/o .gz>
#    Return:
#        dict of original names:filtered names
#        E.g. {'bedgraph/profile.bedGraph.gz': 'wdir/x_profile.bedGraph', 'annotation/genes.gtf.gz': 'wdir/x_genes.gtf'}
#        Output files are *sorted*
#    """
#    nonbam_dict= {}
#    for nonbam in nonbamlist:
#        bname= re.sub('\.gz$', '', os.path.split(nonbam)[1])
#        x_name= os.path.join(tmpdir, 'x_' + bname)
#        pynonbam= pybedtools.BedTool(nonbam)
#        nonbam_x_inbed= pybedtools.BedTool().intersect(a= pynonbam, b= inbed).sort().saveas(x_name)
#        nonbam_dict[nonbam]= x_name
#    return(nonbam_dict)
