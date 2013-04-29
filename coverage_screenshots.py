#!/usr/bin/env python

import argparse
import sys
import pybedtools
import os
import tempfile
import pysam
import subprocess
import shutil
import glob

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

    cat actb.bed
    >chr7	5566757	5566829
    
REQUIREMENTS:
    python 2.7 (2.6 should do, 3.x not sure)
    pysam
    pybedtools
    PyPDF2 (optional: Used to concatenate PDFs in a single one)
    R on $PATH
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--ibam', '-i',
                   required= True,
                   nargs= '+',
                   help='''List of bam files, sorted and indexed to visualize.
Metacharacters are expanded by python (`glob`). E.g. to match all the bam files use
'*.bam'. 
                   ''')

parser.add_argument('--bed', '-b',
                   required= True,
                   help='''Bed file with regions to plot 
                   ''')

parser.add_argument('--fasta', '-f',
                   help='''Fasta file of the reference genome. If given high
resolution plots will show the reference bases at each position.
                   ''')

parser.add_argument('--gtf', '-g',
                   help='''GTF file to fetch annotation from. E.g. gene.gtf in
iGenomes. Can be gzip compressed.
                   ''')

parser.add_argument('--outdir', '-d',
                   required= False,
                   default= None,
                   help='''Output directory for the pdf files. It will be created
if it doesn't exist. Default to current dir.
NB: Not to be confused with --tmpdir where temp file go.
                   ''')

parser.add_argument('--onefile', '-o',
                   required= False,
                   default= None,
                   help='''Concatenate the output PDF files into a single one
passed as argument.
                   ''')

parser.add_argument('--rpm',
                   action= 'store_true',
                   help='''Normalize counts by reads per million using library
sizes. Default is to use raw counts. 
                   ''')

parser.add_argument('--max_nuc', '-m',
                    default= 100,
                    type= int,
                    help='''The maximum width of the region (bp) to plot each nucleotide in
different colour. Default 100 (i.e. regions smaller than 100 bp will be printed with colour
coded nucleotides). 
                   ''')

parser.add_argument('--pheight', '-H',
                    default= 6,
                    type= float,
                    help='''Height of *each* plot in cm. Default 6
                   ''')

parser.add_argument('--pwidth', '-W',
                    default= 24,
                    type= float,
                    help='''Width of the plots in cm. Default 24
                   ''')

parser.add_argument('--psize', '-p',
                    default= 10,
                    type= float,
                    help='''Pointsize for R pdf() function. Sizes
between 9 and 12 should suite most cases. Default 10.
                   ''')

parser.add_argument('--tmpdir', '-t',
                    default= None,
                    help='''Directory where to dump temporary files. If None (default)
python will find one which will be deleted at the end of the execution.
If set it will not be deleted.
                   ''')

# -----------------------------------------------------------------------------
def getFileList(files):
    """Expand the list of files using glob. Return a list of unique files.
    """
    bamlist_dup= []
    for bam in files:
        bamlist_dup.extend(glob.glob(bam))
    bamlist= []
    for bam in bamlist_dup:
        if bam not in bamlist:
            bamlist.append(bam)
    return(bamlist)

def getLibrarySizes(bams):
    """Get the number of reads for each bam file (library sizes)
    bams:
        List of bams for which to get lib sizes
    Returns:
        Dict as {<bam name>:<tot reads>}    
    """
    libsizes= {}
    for bam in bamlist:
        idxstat= pysam.idxstats(bam)
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
    """Compile a string suitable for eval() to execute pysam.mpileup().
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

def pileupBaseCallsToNucs(bases, refbase):
        """Parses the string of read bases from mpileup output to return the count
        of A, C, T, G, N and sum of them. Strandess ignored.
        refbase:
            The reference base for this position. This is the 3rd column in mpileup.
        Return:
            Dictionary with where values are counts and keys cnt_A, cnt_C, cnt_G,
            cnt_T, cnt_N, cnt_Z (this latter being the sum).
        
        See also:
            http://samtools.sourceforge.net/pileup.shtml 
        """
        nuc_counts= {}
        ## Remove the char after '^' since this is mapping quality not base.
        refbase= refbase.upper()
        qscore= [i+1 for i, ltr in enumerate(bases) if ltr == '^']
        calls= [bases[i] for i in range(0, len(bases)) if i not in qscore]
        calls= bases.upper()        
        a= calls.count('A')
        c= calls.count('C')
        g= calls.count('G')
        t= calls.count('T')
        n= calls.count('N')
        r= calls.count('.') + calls.count(',')
        if refbase == 'A':
            a += r
        elif refbase == 'C':
            c += r
        elif refbase == 'G':
            g += r
        elif refbase == 'T':
            t += r
        else:
            n += r
        nuc_counts['cnt_A']= a
        nuc_counts['cnt_C']= c
        nuc_counts['cnt_G']= g
        nuc_counts['cnt_T']= t
        nuc_counts['cnt_N']= n
        nuc_counts['cnt_Z']= a + c + g + t + n
        return(nuc_counts)
        
def parse_pileup(pileup_line, bams):
    """Parse a pileup line (str) typically returned by pysam.mpileup().
    bams:
        List of bam files. Must be in the same order as in mpileup!
    Return:
        Dict with keys: {'chrom': <str>, 'pos': <int>, 'base': <str>,
            <bam.1>: {'depth': int, 'cnt_A': int, 'cnt_C': int, 'cnt_G': int, 'cnt_T': int, 'cnt_N': int},
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
    """Convert the dictionary produced by parse_pileup to a list suitabke to
    be written as bedfile. The bed line as:
    <chrom> <pos-1> <pos> <refbase> <.> <.> <depth.1> <depth.2> ... <depth.n>
    bams:
        list of bam files which are the keys of the dict. MUST be in the same
        order as the list used for parse_pileup()
    """
    bedlist= [pdict['chrom'], pdict['pos']-1, pdict['pos'], pdict['base'], '.', '.']
    for bam in bams:
        bedlist.append(pdict[bam]['depth'])
    for bam in bams:
        bedlist.append(pdict[bam]['cnt_A'])
    for bam in bams:
        bedlist.append(pdict[bam]['cnt_C'])
    for bam in bams:
        bedlist.append(pdict[bam]['cnt_G'])
    for bam in bams:
        bedlist.append(pdict[bam]['cnt_T'])
    for bam in bams:
        bedlist.append(pdict[bam]['cnt_N'])
    for bam in bams:
        bedlist.append(pdict[bam]['cnt_Z'])
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
        List of lists with inner list ['chrom', 'pos', 'base']
    """
    bedregion= pybedtools.BedTool(str(region), from_string= True)
    seq = bedregion.sequence(fi=fasta, tab= True)
    seq= open(seq.seqfn).read().split('\t')
    seq_table= zip([region.chrom] * (region.end - region.start), range(region.start+1, region.end+1), list(seq[1].strip()))
    return(seq_table)
    
def RPlot(**kwargs):
    rout= open(kwargs['rscript'], 'w')
    rplot= """#!/usr/bin/env Rscript
makeTransparent<-function(someColor, alpha=100){
    "Given a colour name (e.g. 'red'), make it transparent.
    someColor:
    Vector of colour names to make transparent e.g. c('red', 'blue')
    alpha:
    Alpha transparency. 100 fully opaque, 0 fully transparent.
    Credit: http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
    "
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

## Graphical parameters
lwd= 4 ## Line width
colList<- list(
    colA= makeTransparent('green', 95),
    colC= makeTransparent('blue', 95),
    colG= makeTransparent('orange', 95),
    colT= makeTransparent('red', 95),
    colZ= 'grey'
)
tgrey<- makeTransparent('blue', 80)

## INPUT
## First read header line only as data to remove paths
header<- read.table('%(mcov)s', header= FALSE, sep= '\t', stringsAsFactors= FALSE, nrows= 1, comment.char= '')
header<- sapply(header[1,], basename)

## Data: Do not read header. It will be added.
mcov<- read.table('%(mcov)s', header= FALSE, sep= '\t', stringsAsFactors= FALSE, skip= 1, comment.char= '')
names(mcov)<- header

## Reference bases
refbases<- read.table('%(refbases)s', header= TRUE, sep= '\t', stringsAsFactors= FALSE, comment.char= '')

## GTF file
do.gtf<- FALSE
if('%(gtf)s' != ''){
    gtf<- read.table('%(gtf)s', header= TRUE, sep= '\t', stringsAsFactors= FALSE, comment.char= '')
    if(nrow(gtf) > 0){
        do.gtf<- TRUE
    }
}

## Column indexes of the counts
count_pos<- list(
    cnt_Z= grep('\\\.cnt_Z$', names(mcov), perl= TRUE), ## Indexes of columns with sum of A+C+G+T
    cnt_A= grep('\\\.cnt_A$', names(mcov), perl= TRUE),
    cnt_C= grep('\\\.cnt_C$', names(mcov), perl= TRUE),
    cnt_G= grep('\\\.cnt_G$', names(mcov), perl= TRUE),
    cnt_T= grep('\\\.cnt_T$', names(mcov), perl= TRUE)
)

## If the range is too wide, use only one colour:
region_size<- max(mcov$end) - min(mcov$start)
if(region_size > %(max_nuc)s) {
    colList<- lapply(colList, function(x) return('grey'))
}

xpos<- rowMeans(mcov[, c('start', 'end')]) + 0.5
nplots<- %(nplots)s
ylim<- c(0, max(mcov[, 4:ncol(mcov)]))

pdf('%(pdffile)s', width= %(pwidth)s/2.54, height= (%(pheight)s*nplots)/2.54, pointsize= %(psize)s)
par(mfrow= c(nplots, 1), las= 1, mar= c(0.5, 4, 0.5, 1), oma= c(3, 1, 3, 0), bty= 'l', mgp= c(3, 0.7, 0))
for(p in seq(1, nplots)){
    libname<- sub('\\\.bam\\\.depth', '', names(mcov)[p+3], perl= TRUE)
    Z<- mcov[, count_pos$cnt_Z[p]]
    A<- mcov[, count_pos$cnt_A[p]]
    C<- mcov[, count_pos$cnt_C[p]] + A
    G<- mcov[, count_pos$cnt_G[p]] + C
    T<- mcov[, count_pos$cnt_T[p]] + G
    plot(xpos, Z, type= 'n', xlab= '', ylab= '', xaxt= 'n', ylim= ylim, lwd= lwd, xlim= c(%(xlim1)s, %(xlim2)s))
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= "grey95", border= 'transparent')
    grid(col= 'darkgrey')
    rect(xleft= mcov$start, ybottom= rep(0, length(xpos)), xright= mcov$end, ytop= Z, col= colList$colZ, border= 'transparent')
    rect(xleft= mcov$start, ybottom= rep(0, length(xpos)), xright= mcov$end, ytop= T, col= colList$colT, border= 'transparent')
    rect(xleft= mcov$start, ybottom= rep(0, length(xpos)), xright= mcov$end, ytop= G, col= colList$colG, border= 'transparent')
    rect(xleft= mcov$start, ybottom= rep(0, length(xpos)), xright= mcov$end, ytop= C, col= colList$colC, border= 'transparent')
    rect(xleft= mcov$start, ybottom= rep(0, length(xpos)), xright= mcov$end, ytop= A, col= colList$colA, border= 'transparent')
    mtext(side= 3, text= libname, adj= 0.02, line= -1, col= tgrey, cex= 0.9)
    if(p == 1){
        if(do.gtf){
            yTop= par('usr')[4]
            y0 <- yTop * 1.05
            arrows(x0= gtf$start, y0= y0, x1= gtf$end, y1= y0,
                col= 'firebrick4', ## ifelse(gtf$type == 'CDS', 'firebrick4', gtf$col),
                lwd= ifelse(gtf$type == 'CDS', 4, 1),
                length= 0.05,
                code= ifelse(gtf$strand == '+', 2, ifelse(gtf$strand == '-', 1, 0)), xpd= NA)
            text(labels= gtf$name, x= gtf$start, y= yTop * 1.1, ifelse(nplots > 3, 0.7, 0.8), col= 'darkgrey', xpd= NA)
        }
    }
}
x<- axis(side= 1, labels= FALSE)
axis(labels= formatC(x, format= 'd', big.mark= ','), side= 1, at= x)
#mtext(text= '%(plotname)s', cex= 1, outer= TRUE, side= 3, xpd= NULL)
mtext(text= '%(plotname)s', cex= 0.95, outer= TRUE, side= 2, las= 0, line= -0.2, adj= 1, col= 'grey50')
mtext(text= '%(ylab)s', cex= 0.95, outer= TRUE, side= 2, las= 0, line= -0.2)
if(nrow(refbases) > 0){
    mtext(at= refbases$pos, side= 1, text= refbases$base, line= 2, cex= ifelse(nplots > 3, 0.66, 0.75), adj= 1, font= 11)
}
dev.off()
""" %kwargs ## %{'mcov': tmp.name, 'plotname': plotname + '.pdf', 'nplots': len(bamlist), 'pheight': args.pheight, 'pwidth':args.pwidth, 'psize': args.psize}
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
    output = PyPDF2.PdfFileWriter()
    for pdf in in_pdf:
        pdfOne = PyPDF2.PdfFileReader(file(pdf, "rb"))
        output.addPage(pdfOne.getPage(0))
    outputStream = file(out_pdf, "wb")
    output.write(outputStream)
    outputStream.close()

def prepare_annotation(gtf_file, bedinterval, outfile):
    """Output an annotation file suitable for R plotting.
    gtf_file:
        Annotation file in gtf format. E.g. gene.gtf for hg19 in iGenomes.
    bedinterval:
        A pybedtools.Interval to intersect to the gtf. Intersected feature will be
        sent to output
    outfile:
        Name for output file
    """
    outf= open(outfile, 'w')
    header= '\t'.join(['chrom', 'start', 'end', 'name', 'type', 'strand', 'col', 'lwd'])
    outf.write(header + '\n')
    gtf= pybedtools.BedTool(gtf_file)
    tmp= tempfile.NamedTemporaryFile(delete= False, suffix= '.bed.txt', prefix= 'gtf_')
    tmp.write(str(bedinterval))
    tmp.close()
    region_gtf= gtf.intersect(pybedtools.BedTool(tmp.name)) ## BedInterval has to be convetred to BedTool
    os.remove(tmp.name)
    feature_graph= {'CDS':{'col': 'red', 'lwd': 4},
                    'exon':{'col': 'blue', 'lwd': 2},
                    'start_codon':{'col': 'blue', 'lwd': 4},
                    'stop_codon':{'col': 'blue', 'lwd': 4}
    }
    for line in region_gtf:
        ftype= line.fields[2]
        pline= [line.chrom, line.start, line.end, line.name, ftype, line.strand]
        pline.extend([feature_graph[ftype]['col'], feature_graph[ftype]['lwd']])
        outf.write('\t'.join([str(x) for x in pline]) + '\n')
    outf.close()
# -----------------------------------------------------------------------------
NWINDS= 1000 ## Number of windows to divide each bed region. Regions smaller than
             ## this will be plotted at base resolution. Larger regions will be
             ## divided in ~1000 equally sized windows and the counts from mpileup
             ## averaged.

def main():
    args = parser.parse_args()
    bamlist= getFileList(args.ibam)
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
    
    ## Make header line for grouped bed files (*.grp.bed.txt)
    ## --------------------------------------------------------------------------
    header= ['chrom', 'start', 'end']
    header.extend([x + '.depth' for x in bamlist])
    header.extend([x + '.cnt_A' for x in bamlist])
    header.extend([x + '.cnt_C' for x in bamlist])
    header.extend([x + '.cnt_G' for x in bamlist])
    header.extend([x + '.cnt_T' for x in bamlist])
    header.extend([x + '.cnt_N' for x in bamlist])
    header.extend([x + '.cnt_Z' for x in bamlist])
    header= '\t'.join(header)
    
    outputPDF= [] ## List of all the pdf files generated. Used only for --onefile
    ## -------------------------------------------------------------------------
    if args.rpm:
        sys.stdout.write('Getting library sizes... ')
        libsizes_dict= getLibrarySizes(bamlist)
        libsizes= [libsizes_dict[x] for x in bamlist]
        print(', '.join([str(x) for x in libsizes]))
    
    inbed= pybedtools.BedTool(args.bed)
    for region in inbed:
        print('Processing: %s' %(str(region).strip()))
        regname= region.chrom + '_' + str(region.start) + '_' + str(region.end)
        ## Reference  FASTA file
        ## File with reference sequence. Open it even if it is going to be header only.
        ## ---------------------------------------------------------------------
        fasta_seq_name= os.path.join(tmpdir, regname + '.seq.txt')
        region_seq= open(fasta_seq_name, 'w')
        region_seq.write('\t'.join(['chrom', 'pos', 'base']) + '\n')
        if ((region.end - region.start) <= 100) and args.fasta:
            fasta_seq= getRefSequence(args.fasta, region)
            for line in fasta_seq:
                region_seq.write('\t'.join([str(x) for x in line]) + '\n')
        region_seq.close()
        ## Get annotation:
        ## --------------------------------------------------------------------
        annot_file= ''
        if args.gtf:
            annot_file= os.path.join(tmpdir, regname + '.annot.txt')
            prepare_annotation(args.gtf, region, annot_file)
        ## Pileup: This is the time consuming part
        ## --------------------------------------------------------------------
        r= '-r ' + region.chrom + ':' + str(region.start) + '-' + str(region.end)
        mpile_cmd= compile_mpileup(bamlist, '-BQ0', '-d1000000', r)
        mpileup_name= os.path.join(tmpdir, regname) + '.mpileup.bed.txt'
        mpileup_bed= open(mpileup_name, 'w')
        for p in eval(mpile_cmd):
            pd= parse_pileup(p, bamlist)
            bedline= pileupToBed(pd, bamlist)
            if args.rpm:
                bedline= bedline[0:6] + rpm(bedline[6:], libsizes)
            mpileup_bed.write('\t'.join([str(x) for x in bedline]) + '\n')
        mpileup_bed.close()
        ## Divide interval in this many regions. No difference if region span < NWINDS
        ## --------------------------------------------------------------------
        w= makeWindows(region, NWINDS) 
        ## Assign to each pileup position its window --------------------------
        mpileup_winds= pybedtools.BedTool(mpileup_name).intersect(w, wb= True)
        ## Aggregate counts in each position by window: -----------------------
        pile_cols= range(7, len(bedline)+1) ## Indexes of columns with counts
        wind_idx= [len(bedline)+1, len(bedline)+2, len(bedline)+3] ## These are the indexes of the columns containing the windows
        mpileup_grp= mpileup_winds.groupby(g= wind_idx, c= pile_cols, o= ['mean'] * len(pile_cols), stream= False)
        mpileup_grp_name= os.path.join(tmpdir, regname) + '.grp.bed.txt'
        mpileup_grp_fout= open(mpileup_grp_name, 'w')
        mpileup_grp_fout.write(header + '\n')
        for line in mpileup_grp:
            mpileup_grp_fout.write(str(line))
        mpileup_grp_fout.close()
        # ----------------------------------------------------------------------
        # Plotting 
        # ----------------------------------------------------------------------
        pdffile= os.path.join(tmpdir, regname + '.pdf')
        outputPDF.append(pdffile)
        rscript= os.path.join(tmpdir, regname + '.R')
        if args.rpm:
            ylab= 'Reads per million'
        else:
            ylab= 'Read count'
        ## Memo: All the args passed to RPlot() become part of the R script.
        rgraph= RPlot(pdffile= pdffile,
              rscript= rscript,
              plotname= regname,
              mcov= mpileup_grp_name,
              refbases= fasta_seq_name,
              nplots= len(bamlist),
              pheight= args.pheight,
              pwidth= args.pwidth,
              psize= args.psize,
              ylab= ylab,
              xlim1= region.start,
              xlim2= region.end,
              gtf= annot_file,
              max_nuc= args.max_nuc)
        if rgraph['stderr'] != '':
            print('\nExpection in executing R script "%s"\n' %(rscript))
            print(rgraph['stdout'])
            sys.exit(rgraph['stderr'])
        if not onefile:
            ## Copy PDFs from temp dir to output dir. Unless you wnated them in onefile
            shutil.copyfile(pdffile, os.path.join(outdir, regname + '.pdf'))
    if onefile:
        catPdf(in_pdf= outputPDF, out_pdf= args.onefile)
    if args.tmpdir is None:
        shutil.rmtree(tmpdir)
if __name__ == '__main__':
    main()
    sys.exit()
