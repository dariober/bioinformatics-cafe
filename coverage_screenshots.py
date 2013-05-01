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

    """, formatter_class= argparse.RawTextHelpFormatter)

# -----------------------------------------------------------------------------
input_args= parser.add_argument_group('Input options', '')

input_args.add_argument('--ibam', '-i',
                   required= False,
                   default= [],
                   nargs= '+',
                   help='''List of bam files, sorted and indexed to visualize.
Metacharacters are expanded by python (`glob`). E.g. to match all the bam files use
'*.bam'. 
                   ''')

input_args.add_argument('--bed', '-b',
                   required= True,
                   help='''Bed file with regions to plot 
                   ''')

input_args.add_argument('--fasta', '-f',
                   help='''Fasta file of the reference genome. If given, high
resolution plots will show the reference bases at each position.
                   ''')

input_args.add_argument('--gtf', '-g',
                    help='''GTF file to fetch annotation from. E.g. gene.gtf in
iGenomes. Can be gzip compressed.
                   ''')

# -----------------------------------------------------------------------------
output_args = parser.add_argument_group('Output options', '')

output_args.add_argument('--outdir', '-d',
                   required= False,
                   default= None,
                   help='''Output directory for the pdf files. It will be created
if it doesn't exist. Default to current dir.
NB: Not to be confused with --tmpdir where temp file go.
                   ''')

output_args.add_argument('--onefile', '-o',
                   required= False,
                   default= None,
                   help='''Concatenate the output PDF files into a single one
passed as argument.
                   ''')

output_args.add_argument('--tmpdir', '-t',
                    default= None,
                    help='''Directory where to dump temporary files. If None (default)
python will find one which will be deleted at the end of the execution.
If set it will not be deleted.
                   ''')

output_args.add_argument('--replot',
                   action= 'store_true',
                   help='''Re-use the output files from a previous execution. Just
redraw the plots using different graphical parameters.
This option allows to reformat the plots without going thorugh the time consuming
steps of generating pileups etc. All but graphical parameters must be the same
as in the execution that generated the original data (obviously --tmpdir must not
be None!).
                   ''')

output_args.add_argument('--rpm',
                   action= 'store_true',
                   help='''Normalize counts by reads per million using library
sizes. Default is to use raw counts. 

                   ''')

output_args.add_argument('--verbose', '-v',
                   action= 'store_true',
                   help='''Print verbose output. Currently this option only adds the
stdout and stderr from R.

                   ''')
# -----------------------------------------------------------------------------
graph_args = parser.add_argument_group('Graphical parameters',
    '''These options passed to R script''')

graph_args.add_argument('--ylim', '-y',
                    default= 'max',
                    type= str,
                    help='''How the maximum value for the y-axis should be set.
Options are: 'max' (default) all y-axes set to the maximum value of all the plots.
'indiv': Scale each plot individually to its maximum
<float>: Set all the plots to this maximum (E.g. 1000).
The lower limit of the y-axis is always 0.
                   ''')

graph_args.add_argument('--col_cov', default= 'grey', help='''Colour for coverage of pooled bases and for N''')
graph_args.add_argument('--col_nuc', nargs= 4, default= '', help='''List of 4 R colours for bars of A, C, G, and T.''')
graph_args.add_argument('--bg', default= 'grey95', help='''Colour for plot background''')
graph_args.add_argument('--nogrid', action= 'store_true', help='''Do not plot grid''')
graph_args.add_argument('--col_text_ann', default= 'black', help='''Colour for annotation text (gene names)''')
graph_args.add_argument('--col_ann', default= 'firebrick4', help='''Colour for annotation bars (exons, CDS etc.)''')
graph_args.add_argument('--col_name', default= '#0000FF50', help='''Colour for text of the name of each plot. Default makeTransparent('blue', 80)''')
graph_args.add_argument('--cex_axis', default= -1, type= float, help='''Character exapansion for the axis annotation (cex.axis in R).
Use negative value to set default. ''')
graph_args.add_argument('--cex_range', default= -1, type= float, help='''Character exapansion for the text range of plot''')
graph_args.add_argument('--cex_seq', default= -1, type= float, help='''Character exapansion for the nucleotide sequence''')
graph_args.add_argument('--line_range', default= 2, type= float, help='''Distance of range bar from x-axis. In R's line units''')
graph_args.add_argument('--line_seq', default= 3.5, type= float, help='''Distance of nucleotide sequence bar from x-axis. In R's line units''')
graph_args.add_argument('--col_seq', default= 'black', help='''Colour for the nucleotide sequence''')
graph_args.add_argument('--oma', default= [5, 1.1, 3, 1.1], nargs= 4, type= float, help='''List of 4 floats giving the outer margins of the plot.
Default 4 1.1 3 1.1''')
graph_args.add_argument('--mar', default= [0.5, 4, 0.5, 1], nargs= 4, type= float, help='''List of 4 floats giving the margins of each plot.
Default 0.5 4 0.5 1''')

graph_args.add_argument('--max_nuc', '-m',
                    default= 100,
                    type= int,
                    help='''The maximum width of the region (bp) to plot each nucleotide in
different colour. Default 100 (i.e. regions smaller than 100 bp will be printed with colour
coded nucleotides). 
                   ''')

graph_args.add_argument('--max_fa', '-M',
                    default= 100,
                    type= int,
                    help='''The maximum width of the region (bp) to plot each nucleotide at the
bottome of the lower plot. Default 100 (i.e. regions smaller than 100 bp will have the refence
sequence shown). Ignored if --fasta is not given 
                   ''')

graph_args.add_argument('--pheight', '-H',
                    default= -1,
                    type= float,
                    help='''Height of *each* plot in cm. Default 6
                   ''')

graph_args.add_argument('--pwidth', '-W',
                    default= 15,
                    type= float,
                    help='''Width of the plots in cm. Default 24
                   ''')

graph_args.add_argument('--psize', '-p',
                    default= 10,
                    type= float,
                    help='''Pointsize for R pdf() function. Sizes
between 9 and 12 should suite most cases. Default 10.
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
    for bam in bams:
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
        for x in bases:
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
    """Parse a pileup line (str) typically returned by pysam.mpileup().
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
    """Convert the dictionary produced by parse_pileup to a list suitabke to
    be written as bedfile. The bed line as:
    <chrom> <pos-1> <pos> <refbase> <.> <.> <depth.1> <depth.2> ... <depth.n>
    bams:
        list of bam files which are the keys of the dict. MUST be in the same
        order as the list used for parse_pileup()
    """
    bedlist= [pdict['chrom'], pdict['pos']-1, pdict['pos'], pdict['base'], '.', '.']
    for c in ['depth', 'A', 'C', 'G', 'T', 'N', 'Z']:
        for bam in bams:
            bedlist.append(pdict[bam][c])
#    for bam in bams:
#        bedlist.append(pdict[bam]['depth'])
#    for bam in bams:
#        bedlist.append(pdict[bam]['A'])
#    for bam in bams:
#        bedlist.append(pdict[bam]['C'])
#    for bam in bams:
#        bedlist.append(pdict[bam]['G'])
#    for bam in bams:
#        bedlist.append(pdict[bam]['T'])
#    for bam in bams:
#        bedlist.append(pdict[bam]['N'])
#    for bam in bams:
#        bedlist.append(pdict[bam]['Z'])
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
# ------------------------------------------------------------------------------
# Intial settings
# ------------------------------------------------------------------------------
lwd= 4
cex.axis<- ifelse(%(cex_axis)s < 0, par('cex.axis'), %(cex_axis)s)

col_nuc<- c(%(col_nuc)s)
colList<- list(
    colA= ifelse(is.null(col_nuc), makeTransparent('green', 95), col_nuc[1]),
    colC= ifelse(is.null(col_nuc), makeTransparent('blue', 95), col_nuc[2]),
    colG= ifelse(is.null(col_nuc), makeTransparent('orange', 95), col_nuc[3]),
    colT= ifelse(is.null(col_nuc), makeTransparent('red', 95), col_nuc[4]),
    colZ= '%(col_cov)s'
)

# ------------------------------------------------------------------------------
# INPUT
# ------------------------------------------------------------------------------
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
    gtf<- gtf[which(gtf$type %%in%% c('start_codon', 'stop_codon') == FALSE),] ## Do not annotate start and stop codon
    if(nrow(gtf) > 0){
        do.gtf<- TRUE
    }
}

# ------------------------------------------------------------------------------
# Column indexes of the counts
count_pos<- list(
    Z= grep('\\\.Z$', names(mcov), perl= TRUE), ## Indexes of columns with sum of A+C+G+T+N
    A= grep('\\\.A$', names(mcov), perl= TRUE),
    C= grep('\\\.C$', names(mcov), perl= TRUE),
    G= grep('\\\.G$', names(mcov), perl= TRUE),
    T= grep('\\\.T$', names(mcov), perl= TRUE)
)
## A check: all the counts above have the same length:
if(length(unique((sapply(count_pos, length)))) != 1){
    stop('An error occured while processing the data.')
}

## If the range is too wide, use only one colour:
region_size<- max(mcov$end) - min(mcov$start)
if(region_size > %(max_nuc)s) {
    colList<- lapply(colList, function(x) return(colList$colZ))
}

xpos<- rowMeans(mcov[, c('start', 'end')]) + 0.5
nplots<- length(count_pos[['Z']])

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
pwidth<- %(pwidth)s
pheight<- %(pheight)s
if(pheight <= 0){
    #Get sensible default values for height
    pheight<- (pwidth * 0.35) + ((pwidth/10) / nplots)
}
pdf('%(pdffile)s', width= pwidth/2.54, height= (pheight * nplots)/2.54, pointsize= %(psize)s)
par(mfrow= c(nplots, 1), las= 1, mar= c(%(mar)s), oma= c(%(oma)s), bty= 'l', mgp= c(3, 0.7, 0))
for(p in seq(1, nplots)){
    ## For each library:
    ## -----------------
    libname<- sub('\\\.bam\\\.depth', '', names(mcov)[p+3], perl= TRUE)
    Z<- mcov[, count_pos$Z[p]]
    A<- mcov[, count_pos$A[p]]
    C<- mcov[, count_pos$C[p]] + A
    G<- mcov[, count_pos$G[p]] + C
    T<- mcov[, count_pos$T[p]] + G

    ## Set maximum for y-axt
    ## ---------------------
    if('%(ylim)s' == 'max'){
        ylim<- max(mcov[, 5:ncol(mcov)])
    } else if('%(ylim)s' == 'indiv'){
        ylim<- max(Z)
    } else {
        ylim<- as.numeric('%(ylim)s')
    }
    
    plot(xpos, Z, type= 'n', xlab= '', ylab= '', xaxt= 'n', ylim= c(0, ylim), lwd= lwd, xlim= c(%(xlim1)s, %(xlim2)s), cex.axis= cex.axis)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= "%(bg)s", border= 'transparent')
    if('%(nogrid)s' == 'False'){
        grid(col= 'darkgrey')
    }
    rect(xleft= mcov$start, ybottom= rep(0, length(xpos)), xright= mcov$end, ytop= Z, col= colList$colZ, border= 'transparent')
    rect(xleft= mcov$start, ybottom= rep(0, length(xpos)), xright= mcov$end, ytop= T, col= colList$colT, border= 'transparent')
    rect(xleft= mcov$start, ybottom= rep(0, length(xpos)), xright= mcov$end, ytop= G, col= colList$colG, border= 'transparent')
    rect(xleft= mcov$start, ybottom= rep(0, length(xpos)), xright= mcov$end, ytop= C, col= colList$colC, border= 'transparent')
    rect(xleft= mcov$start, ybottom= rep(0, length(xpos)), xright= mcov$end, ytop= A, col= colList$colA, border= 'transparent')
    mtext(side= 3, text= libname, adj= 0.02, line= -1, col= '%(col_name)s', cex= 0.9)
    if(p == 1){
makeTransparent('blue', 80)
        ## Plotting of annotation
        ## ----------------------
        if(do.gtf){
            yTop= par('usr')[4]
            y0 <- yTop * 1.05
            thick_offset<- y0 * 0.025 ## Amount to add and subtract to y0 to draw thich boxes (CDSs)
            thin_offset<- y0 * 0.015 ## Amount to add and subtract to y0 to draw thich boxes (CDSs)
            rect(xleft= gtf$start,
                        ybottom= y0 - ifelse(gtf$type == 'CDS', thick_offset, thin_offset),
                        xright= gtf$end,
                        ytop= y0 + ifelse(gtf$type == 'CDS', thick_offset, thin_offset),
                        col= '%(col_ann)s', xpd= NA, border= 'transparent')
            ## Plotting of gene names with strand
            gene_strand<- paste(gtf$name, ifelse(gtf$strand %%in%% c('+', '-'), gtf$strand, ''))
            text(labels= gene_strand, x= rowMeans(gtf[,c('start', 'end')]), y= y0 + thick_offset + thin_offset, cex= cex.axis * 0.9, col= '%(col_text_ann)s', xpd= NA, adj= c(0.5,0)) ## 
        }
    }
}

x<- axis(side= 1, labels= FALSE)
axis(labels= formatC(x, format= 'd', big.mark= ','), side= 1, at= x, cex.axis= cex.axis)
mtext(text= '%(plotname)s', cex= 0.95, outer= TRUE, side= 4, las= 0, line= 0, col= 'grey50')
mtext(text= '%(ylab)s', cex= 0.95, outer= TRUE, side= 2, las= 0, line= -0.2)
if(nrow(refbases) > 0){
    cex_seq<- ifelse(%(cex_seq)s <= 0, ifelse(nplots > 3, 0.66, 0.75), %(cex_seq)s)
    mtext(at= refbases$pos,
        side= 1,
        text= refbases$base,
        line= %(line_seq)s,
        cex= cex_seq,
        col= '%(col_seq)s',
        adj= 1,
        family= 'mono',
        font= 1)
}
## Text for range
wcex_range<- ifelse(%(cex_range)s <= 0, ifelse(nplots > 3, 0.66, 0.7), %(cex_range)s)
xrange<- x[length(x)] - x[1]
mtext(text= '|', at= c(x[1], x[length(x)]), line= %(line_range)s, side= 1, xpd= NA, cex= wcex_range)
mtext(text= formatC(paste(xrange, 'bp'), format= 'd', big.mark= ','), line= %(line_range)s, side= 1, xpd= NA, cex= wcex_range)
dev.off()
""" %kwargs ## %{'mcov': tmp.name, 'plotname': plotname + '.pdf', 'pheight': args.pheight, 'pwidth':args.pwidth, 'psize': args.psize}
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

def make_dummy_mpileup(chrom, start, end, nbams):
    """Create an empty line from mpileup to be used for regions w/o any reads in
    any library. mpileup skips such regions altogheter and they wouldn't b plotted
    otherwise.
    Output will look like below with columns after the 6th with 0:
    ['lambda_gi9626243', 9, 10, 'N', '.', '.', 327, 1, 0, 0, 326, 0, 327]
    Each library occupies 7 columns
    chrom, start, end:
        Chrom and position to fill with zeros
    nbams:
        Number of bam files that would be present
    """
    zeros= [0] * 7 * nbams
    bedline= [chrom, start, end, 'N', '.', '.'] + zeros
    return(bedline)
    
# -----------------------------------------------------------------------------
NWINDS= 1000 ## Number of windows to divide each bed region. Regions smaller than
             ## this will be plotted at base resolution. Larger regions will be
             ## divided in ~1000 equally sized windows and the counts from mpileup
             ## averaged.

def main():
    args = parser.parse_args()
    if args.replot and args.tmpdir is None:
        sys.exit('\nCannot replot without a working (--tmpdir) directory!\n')
    bamlist= getFileList(args.ibam)
    if not args.replot:
        print('\nFiles to analyze (%s found):\n%s\n' %(len(bamlist), ', '.join(bamlist)))
    if len(bamlist) == 0 and not args.replot:
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
    
    ## Make header line for grouped bed files (*.grp.bed.txt)
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
    
    outputPDF= [] ## List of all the pdf files generated. Used only for --onefile
        
    ## -------------------------------------------------------------------------
    if args.rpm and not args.replot:
        sys.stdout.write('Getting library sizes... ')
        libsizes_dict= getLibrarySizes(bamlist)
        libsizes= [libsizes_dict[x] for x in bamlist]
        print(', '.join([str(x) for x in libsizes]))
    
    inbed= pybedtools.BedTool(args.bed)
    for region in inbed:
        print('Processing: %s' %(str(region).strip()))
        regname= '_'.join([str(x) for x in [region.chrom, region.start, region.end]])
        if region.name != '':
            regname = regname + '_' + region.name 

        ## Prepare output file names
        ## --------------------------------------------------------------------
        fasta_seq_name= os.path.join(tmpdir, regname + '.seq.txt')
        if args.gtf:
            annot_file= os.path.join(tmpdir, regname + '.annot.txt')
        else:
            annot_file= ''
        mpileup_name= os.path.join(tmpdir, regname) + '.mpileup.bed.txt'
        mpileup_grp_name= os.path.join(tmpdir, regname) + '.grp.bed.txt'
        pdffile= os.path.join(tmpdir, regname + '.pdf')
        rscript= os.path.join(tmpdir, regname + '.R')
        if not args.replot:
            ## Reference  FASTA file
            ## File with reference sequence. Open it even if it is going to be header only.
            ## ---------------------------------------------------------------------
            region_seq= open(fasta_seq_name, 'w')
            region_seq.write('\t'.join(['chrom', 'pos', 'base']) + '\n')
            if ((region.end - region.start) <= args.max_fa) and args.fasta:
                fasta_seq= getRefSequence(args.fasta, region)
                for line in fasta_seq:
                    region_seq.write('\t'.join([str(x) for x in line]) + '\n')
            region_seq.close()
            ## Get annotation:
            ## --------------------------------------------------------------------
            #annot_file= ''
            if args.gtf:
                prepare_annotation(args.gtf, region, annot_file)
            ## Pileup: This is the time consuming part
            ## --------------------------------------------------------------------
            r= '-r ' + region.chrom + ':' + str(region.start) + '-' + str(region.end)
            mpile_cmd= compile_mpileup(bamlist, '-BQ0', '-d1000000', r)
            mpileup_bed= open(mpileup_name, 'w')
            for p in eval(mpile_cmd):
                pd= parse_pileup(p, bamlist)
                bedline= pileupToBed(pd, bamlist)
                if args.rpm:
                    bedline= bedline[0:6] + rpm(bedline[6:], libsizes)
                mpileup_bed.write('\t'.join([str(x) for x in bedline]) + '\n')
            mpileup_bed.close()
            if os.stat(mpileup_name).st_size == 0:
                mpileup_bed= open(mpileup_name, 'w')
                bedline= make_dummy_mpileup(region.chrom, region.start, region.start + 1, len(bamlist))
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
            mpileup_grp_fout= open(mpileup_grp_name, 'w')
            mpileup_grp_fout.write(header + '\n')
            for line in mpileup_grp:
                mpileup_grp_fout.write(str(line))
            mpileup_grp_fout.close()
        # ----------------------------------------------------------------------
        # Plotting 
        # ----------------------------------------------------------------------
        outputPDF.append(pdffile)
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
              pheight= args.pheight,
              pwidth= args.pwidth,
              psize= args.psize,
              ylab= ylab,
              xlim1= region.start,
              xlim2= region.end,
              gtf= annot_file,
              max_nuc= args.max_nuc,
              ylim= args.ylim,
              cex_axis= args.cex_axis,
              col_cov= args.col_cov,
              col_nuc= quoteStringList(args.col_nuc),
              bg= args.bg,
              nogrid= args.nogrid,
              col_text_ann= args.col_text_ann,
              col_ann= args.col_ann,
              col_name= args.col_name,
              cex_range= args.cex_range,
              cex_seq= args.cex_seq,
              line_range= args.line_range,
              line_seq= args.line_seq,
              col_seq= args.col_seq,
              oma= ', '.join([str(x) for x in args.oma]),
              mar= ', '.join([str(x) for x in args.mar]))
        
        
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
