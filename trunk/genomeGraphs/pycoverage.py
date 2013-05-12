import sys
import pybedtools
import os
import tempfile
import subprocess
import glob
import pycoverage
import gzip
import re

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

def getLibrarySizes(bams, samtools_path= ''):
    """Get the number of reads for each bam file (library sizes)
    bams:
        List of bams for which to get lib sizes
    Returns:
        Dict as {<bam name>:<tot reads>}    
    """
    samtools_idx= os.path.join(samtools_path, 'samtools idxstats')
    libsizes= {}
    for bam in bams:
        cmd= samtools_idx + ' ' + bam
        proc= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
        idxstat, idxerr= proc.communicate()
        idxstat= idxstat.strip().split('\n')
        idxstat= [x.split('\t') for x in idxstat]
        libsize= sum([int(x[2]) for x in idxstat])
        libsizes[bam]= libsize
    return(libsizes)

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

def normMultiCovLine(line):
    """line is a line of output from multi_bam_coverage. Divide each count by
    the interval size to normalize it.
    ['chr1', '99', '100', '1090', '182']
    Return:
        List of strings with the same elements as input with counts normalized
    """
    counts= [float(x) for x in line.fields[3:]]
    intsize= line.end - line.start
    normcounts= line.fields[0:3] + [str(x/intsize) for x in counts]    
    return(normcounts)

def mergePDF(filenames, output_filename):
    output = PyPDF2.PdfFileWriter()
    for filename in filenames:
        input = PyPDF2.PdfFileReader(file(filename, "rb"))
        for page in input.pages:
            output.addPage(page)
    outputstream = file(output_filename, "wb")
    output.write(outputstream)
    outputstream.close()

def compile_mpileup(bams, *args):
    """DEPRECATED: Compile a string suitable for eval() to execute pysam.mpileup().
    Note: syntax `pysam.mpileup(['bam1', 'bam2']) is not supported!`
    bams:
        List of bam names
    *args:
        Further arguments passed to mpileup. E.g '-BQ 0', '-d 10000000'
    Return:
        String to be passed to eval()
    """
    mpile= 'pysam.mpileup('
    for bam in list(args) + bams:
        bamarg= '"' + bam + '", '
        mpile += bamarg
    mpile= mpile[:-2] ## Remove last command & space
    mpile += ')' 
    return(mpile)

def mpileup_cmd(bamlist, region, fasta= None, mpileup= 'samtools mpileup'):
    """Compile a command string to execute samtools mpileup
    bamlist:
        List of input bams
    region:
        Pybedtools interval to get chrom start, end position
    f:
        FASTA file to get sequence from. If file is not indexed.
    mpileup:
        String with the full path to mpileup or just samtools mpileup
        if the rogram is opn path
    Return:
        String
    """
    r= '-r ' + region.chrom + ':' + str(region.start) + '-' + str(region.end)
    if fasta:
        f= '-f %s' %(fasta)
    else:
        f= ''        
    cmd= '%(mpileup)s %(f)s -BQ0 -d10000000 %(r)s %(bamlist)s' %{'mpileup': mpileup, 'f': f, 'r': r, 'bamlist': ' '.join(bamlist)}
    return(cmd)
                
def pileupBaseCallsToNucs(bases, refbase):
        """Parses the string of read bases from mpileup output to return the count
        of A, C, T, G, N and sum of them. Strandess ignored.
        refbase:
            The reference base for this position. This is the 3rd column in mpileup.
        Return:
            Dictionary with where values are counts and keys A, C, G,
            T, N, Z (this latter being the sum).
        
        See also:
            http://samtools.sourceforge.net/pileup.shtml 
        """
        nuc_counts= {}
        ## Remove the char after '^' since this is mapping quality not base.
        refbase= refbase.upper()
        callDict= {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, '.': 0, ',': 0}
        keys= tuple(callDict.keys())
        skip= False
        for x in bases.upper():
            if x  == '^':
                skip= True
            elif skip:
                skip= False
            elif x in keys:
                callDict[x] += 1
            else:
                pass
        callDict[refbase] += (callDict['.'] + callDict[','])
        callDict['Z']= sum((callDict[x] for x in ('A', 'C', 'G', 'T', 'N')))
        return(callDict)
        
def parse_pileup(pileup_line, bams):
    """Parse a pileup line (str) typically returned by pysam.mpileup() or
    mpileup via subprocess.
    bams:
        List of bam files. Must be in the same order as in mpileup!
    Return:
        Dict with keys: {'chrom': <str>, 'pos': <int>, 'base': <str>,
            <bam.1>: {'depth': int, 'A': int, 'C': int, 'G': int, 'T': int, 'N': int},
            <bam.2>: {...}, ...}
    """
    pdict= {}
    plist= pileup_line.split('\t')
    pdict['chrom']= plist[0]
    pdict['pos']= int(plist[1])
    pdict['base']= plist[2]
    N= 3
    for bam in bams:
        pdict[bam]= {'depth': int(plist[N])}
        N += 1
        baseDict= pileupBaseCallsToNucs(plist[N], pdict['base'])
        pdict[bam].update(baseDict) ## Concatenate nuc counts to existing dict. 
        N += 2
    return(pdict)

def pileupToBed(pdict, bams):
    """Convert the dictionary produced by parse_pileup to a list suitable to
    be written as bedfile. The bed line as:
    <chrom> <pos-1> <pos> <depth.1> <depth.2> ... <depth.n>
    bams:
        list of bam files which are the keys of the dict. MUST be in the same
        order as the list used for parse_pileup()
    """
    bedlist= [pdict['chrom'], pdict['pos']-1, pdict['pos']]
    for c in ['depth', 'A', 'C', 'G', 'T', 'N', 'Z']:
        for bam in bams:
            bedlist.append(pdict[bam][c])
    return(bedlist)

def rpm(raw_counts, libsize):
    """Normalize counts by dividing libsize and x1000000 (Reads Per Million)
    raw_counts:
        List of ints of raw counts to normalize
    libsize:
        List of library sizes to divide each count. Recycled
    Return:
        List of floats.
    Example:
    raw_counts= [10,     500,     100,   500]
    libsize=    [10000,  50000]
    rpm=        [1000.0, 10000.0, 10000.0, 10000.0]
    """
    expLibList= libsize * (len(raw_counts) / len(libsize))
    if len(expLibList) != len(raw_counts):
        return(False)
    rpmList= []
    for r,s in zip(raw_counts, expLibList):
        rpmList.append((float(r)/s)*1000000)
    return(rpmList)

def getRefSequence(fasta, region):
    """Read the fasta file and extract the region given in bedtool interval region.
    Return:
        List of tuples with inner tuple ['chrom', 'start', 'end', 'base']
    NB: You need to reduce the start by 1 because fastaFromBed seems to be 1-based.
    """
    region_str= '\t'.join([region.chrom, str(region.start-1), str(region.end)])
    bedregion= pybedtools.BedTool(str(region_str), from_string= True)
    seq = bedregion.sequence(fi=fasta, tab= True)
    seq= open(seq.seqfn).read().split('\t')
    seq_table= zip(
        [region.chrom] * (region.end - (region.start-1)), ## Column of chrom
        range(region.start-1, region.end), ## Column of start pos
        range(region.start, region.end+1), ## Col of end pos
        list(seq[1].strip()) ## Col of bases
        )
    return(seq_table)

def quoteStringList(x):
    '''Join the list of string x in a string where each element is double quoted and
    separted by comma & space.
    Useful to convert python's list of strings to R character vectors
    E.g. x= ['blue', 'black']
    quoteStringList(x) >>> '"blue", "black"'
    NB:
        Numbers are converted to strings and quoted as well!
    NB2:
        This function is not meant to cope well with strings containing double quotes,
        weird metachars etc.
    '''
    s= ''
    for y in x:
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
    rplot= rtemplate %kwargs ## %{'mcov': tmp.name, 'plotname': plotname + '.pdf', 'pheight': args.pheight, 'pwidth':args.pwidth, 'psize': args.psize}
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

def prepare_nonbam_file(infile_name, outfile_handle, region):
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
    Return:
        Number of intervals that overlap `region` (n lines)
    """
    infile= pybedtools.BedTool(infile_name)
    region_x_infile= infile.all_hits(region)
    isGTF= False
    nlines= 0
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
            outline= [line.chrom, line.start - 1, line.end, infile_name, 'NA', 'NA', 'NA', 'NA', 'NA', line.fields[2], name, strand] ## NA for ACTGZ
        elif infile_name.lower().endswith('.bedgraph') or infile_name.lower().endswith('.bedgraph.gz'):
            outline= [line.chrom, line.start, line.end, infile_name, '0', '0', '0', '0', line.name, 'coverage', 'NA', strand]
        else:
            outline= [line.chrom, line.start, line.end, infile_name, 'NA', 'NA', 'NA', 'NA', 'NA', 'generic', line.name, strand]
        outfile_handle.write('\t'.join([str(x) for x in outline]) + '\n')
        nlines += 1
    return(nlines)

def make_dummy_mpileup(chrom, start, end, nbams):
    """Create an empty line from mpileup to be used for regions w/o any reads in
    any library. mpileup skips such regions altogheter and they wouldn't b plotted
    otherwise.
    Output will look like below with columns after the 6th with 0:
    ['lambda_gi9626243', 9, 10, 327, 1, 0, 0, 326, 0, 327]
    Each library occupies 7 columns: depth, A, C, T, G, N, Z.
    chrom, start, end:
        Chrom and position to fill with zeros
    nbams:
        Number of bam files that would be present
    """
    zeros= [0] * 7 * nbams
    bedline= [chrom, start, end] + zeros
    return(bedline)

def prepare_reference_fasta(fasta_seq_name, maxres, region, fasta):
    """Output a reference file with the base at each postion in
    the region interval
    fasta_seq_name:
        Output name for reference file. Will have format <chrom> <pos> <refbase>
        included header.
    maxres:
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
    if ((region.end - region.start) <= maxres) and fasta:
        fasta_seq= getRefSequence(fasta, region)
        for line in fasta_seq:
            region_seq.write('\t'.join([str(x) for x in line]) + '\n')
    region_seq.close()
    return(True)

def bamlist_to_mpileup(mpileup_name, mpileup_grp_name, bamlist, region, nwinds, fasta, RPM, regionWindows, samtools, groupFun= 'mean'):
    """Output mpileup and grouped mpileup files for list of bam files
    mpileup_name, mpileup_grp_name:
        Name for output mpileup and grouped mpileup file
    bamlist:
        List of bam files
    region:
        pybedtools interval where mpileup should be produced (passed to
        -r options of mpileup)
    nwinds:
        Maximum number of positions before windowing kinks in (args.nwinds)
    fasta:
        fasta file for mpileup reference
    RPM:
        True/False for whether mpileup counts should be normlaized to RPM
    regionWindows:
        bed interval divided into regions by (output of) pycoverage.makeWindows
    samtools:
        Path to samtools (just the path, e.g. /home/myself/bin)
    groupFun:
        Apply this function to group-by windows. This opt passed to bedtools
        groupby. Check there for valid options
    Returns:
        True on success. Side effect is to produce *.mpileup.bed.txt, *.grp.bed.txt
    """
    ## Make header line for grouped bed files (*.grp.bed.txt) from mpileup
    ## --------------------------------------------------------------------------
    header= ['chrom', 'start', 'end']
    header.extend([x + '.depth' for x in bamlist])
    header.extend([x + '.A' for x in bamlist])
    header.extend([x + '.C' for x in bamlist])
    header.extend([x + '.G' for x in bamlist])
    header.extend([x + '.T' for x in bamlist])
    header.extend([x + '.N' for x in bamlist])
    header.extend([x + '.Z' for x in bamlist])
    header= '\t'.join(header)

    mpileup_bed= open(mpileup_name, 'w')
    cmd= mpileup_cmd(bamlist= bamlist, region= region, fasta= fasta, mpileup= os.path.join(samtools, 'samtools mpileup'))
    proc= subprocess.Popen(cmd, shell= True, stdout=subprocess.PIPE, stderr= subprocess.PIPE)
    if RPM:
        libsizes= getLibrarySizes(bamlist, samtools_path= samtools)
        libsizes= [libsizes[x] for x in bamlist]
    nlines= 0
    while True:
        ## Use this while loop to avoid reading in memory all the output of mpileup.
        ## See also http://stackoverflow.com/questions/2804543/read-subprocess-stdout-line-by-line
        line= proc.stdout.readline()
        sys.stdout.flush()
        if not line:
            break
        pd= parse_pileup(line, bamlist)
        bedline= pileupToBed(pd, bamlist)
        cnt_indx= 3 ## Column index in bedline where counts start (4th columns)
        if RPM:
            bedline= bedline[0:cnt_indx] + rpm(bedline[cnt_indx:], libsizes)
        mpileup_bed.write('\t'.join([str(x) for x in bedline]) + '\n')
        nlines += 1
    mpileup_bed.close()
    if os.stat(mpileup_name).st_size == 0:
        mpileup_bed= open(mpileup_name, 'w')
        bedline= make_dummy_mpileup(region.chrom, region.start, region.start + 1, len(bamlist))
        mpileup_bed.write('\t'.join([str(x) for x in bedline]) + '\n')
        mpileup_bed.close()
    mpileup_grp_fout= open(mpileup_grp_name, 'w')
    if nlines > nwinds:
        """Divide interval in nwinds regions if the number of positions to plot is >nwinds
        """
        mpileup_winds= pybedtools.BedTool(mpileup_name).intersect(regionWindows, wb= True)  ## Assign to each pileup position its window
        pile_cols= range(cnt_indx+1, len(bedline)+1) ## Indexes of columns with counts 1-BASED!
        wind_idx= [len(bedline)+1, len(bedline)+2, len(bedline)+3] ## These are the indexes of the columns containing the windows
        mpileup_grp= mpileup_winds.groupby(g= wind_idx, c= pile_cols, o= [groupFun] * len(pile_cols), stream= False) ## Aggregate counts in each position by window
        mpileup_grp_fout.write(header + '\n')
        for line in mpileup_grp:
            mpileup_grp_fout.write(str(line))
    else:
        """If all the positions are to be plotted (nlines < nwinds), copy the output of
        mpileup with the header line.
        """
        mpileup_bed= open(mpileup_name)
        mpileup_grp_fout.write(header + '\n')        
        for line in mpileup_bed:
            mpileup_grp_fout.write(str(line))
    mpileup_grp_fout.close()        
    return(True)
    
def compressBedGraph(regionWindows, bedgraph_name, use_file_name, bedgraph_grp_fh, col_idx= 9, groupFun= 'mean'):
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
        outline= [line.chrom, str(line.start), str(line.end), use_file_name, '0', '0', '0', '0', line.name, 'coverage', 'NA', '.']
        bedgraph_grp_fh.write('\t'.join(outline) + '\n')

