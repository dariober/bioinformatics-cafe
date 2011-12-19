#!/usr/bin/python   

"""

Rcode to produce boxplot of mapq scores from bamqc output
--------------------------------------------------------------------------------
indx<- grep('mapq_quant.', names(bamqc_out))
bx<- t(bamqc_out[,indx])
bxp(list(stats= bx, n= NA, conf= NA, group= NA, names= NA))

"""

import pysam
import sys
import os
import argparse
import subprocess

bamfiles= sys.argv[1:]

LIMITS_MAPQ= [0, 5, 10, 15, 20, 25, 30, 35, 40]
LIMITS_NM= [0, 1, 2, 3, 4, 5, 6]
QUANTILES= [0.05, 0.25, 0.5, 0.75, 0.95]

parser = argparse.ArgumentParser(description= """

DESCRIPTION:
                                 
    Perform quality control on one or more bam files.
    Use '-i -' to read the list of input files from standard input.
    
    bamqc.py produces summary statistics related to mapq scores and edit
    distance between read and reference.

OUTPUT:

    The output is:
    - A table with one row per file (output).
    - A png file with plotted quality metrics (output.png)
    - R script with input data and code to produce the above image (output.R)
    - The output outomatically produced by R (bamqc.Rout) 

    Columns in output:
    filename:
        Filename with path stripped.
    len_median, len_sd:
        Read lenght: median and std. dev.
    n:
        Total number of reads
    aln:
        Aligned reads (flag doeas not contain 4)
    perc_aln:
        Percentage read aliged
    mapq.0, mapq.5, ..., mapq.40:
        Number of reads with mapq score 0 or above, 5 or above..., 40 or above     
    mapq_255:
        No. reads with mapq not available (mapq= 255)
    nm.0, nm.1, ..., nm.6:
        Number of reads aligned with edit distance (NM tag) 0, with edit distance 1, 2, 3, 4, 5, 6 or above.
    nm_na:
        Number of reads without NM tag
    mapq_quant.0.05, ..., mapq_quant.0.95:
        Quantiles for mapq score.

EXAMPLES:

    bamqc.py -i myfile1.bam myfile2.bam -o myqc.txt
    
    ## Read files from stdin. Send output to stdout:
    ls *.bam | bamqc.py -i -

TODO:
    - QC specific for PE reads.

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-i', '--input', # nargs='?', default=sys.stdout,
                    type= str,
                    nargs= '+',
                    required= True,
                    help=""" List of input bam files. Use '-i -' to get the list from stdin
                    
                    """)

parser.add_argument('-o', '--output', # nargs='?', default=sys.stdout,
                    type= argparse.FileType('w'),
                    default= sys.stdout,
                    metavar= 'FILE',
                    help="""Output file. Default is stdout
                    
                    """)

parser.add_argument('--nograph', # nargs='?', default=sys.stdout,
                    action= 'store_true',
                    help="""Do not produce R graphics
                    
                    """)
args = parser.parse_args()
outfilename = args.output.name

if args.input == ['-']:
    args.input= sys.stdin.readlines()
    args.input= [x.strip() for x in args.input]


# ------------------------------------------------------------------------------
#                             Classes & Functions
# ------------------------------------------------------------------------------

## Order in which the attributes in class File_Stats should be retured:
header=  ['filename', 'median_length', 'len_sd', 'nreads_tot', 'nreads_aln', 'perc_aln', 'mapq', 'nreads_mapq_255', 'nreads_nm', 'nreads_nm_na', 'mapq_quantiles']

## Column names to output (must reflect the order in 'header' above)
colnames=  ['filename', 'len_median', 'len_sd', 'n', 'aln', 'perc_aln'] +  \
           ['mapq.'+str(x) for x in LIMITS_MAPQ] + \
           ['mapq_255'] + \
           ['nm.'+str(x) for x in LIMITS_NM] + \
           ['nm_na'] + \
           ['mapq_quant.'+str(x) for x in QUANTILES]

class File_Stats:
    def __init__(self):
        """ Attributes here will be columns in output table """
        self.filename= ''
        self.median_length= 'NA'
        self.len_sd= 'NA'
        self.nreads_tot= 0
        self.nreads_aln= 0
        self.perc_aln= 0
        self.nreads_mapq= 0
        self.nreads_mapq_255= 0
        self.nreads_nm= {}
        self.nreads_nm_na= 0
        self.mapq= {}
        self.mapq_quantiles= []

def cumdist(mapq_dict, limits= [0, 5, 10, 15, 20, 25, 30, 35, 40]):
    """
    mapq_dict= {0:10, 4:5, 10:5, 11:1, 15:5, 16:1, 20:3}
    cumdist(mapq_dict)
    >>> [[0, 30], [5, 15], [10, 15], [15, 9], [20, 3], [25, 0], [30, 0], [35, 0], [40, 0]]
    """
    qhist={}
    for q in limits:
        qhist[q]= 0
    for mapq in mapq_dict.keys():
        for histv in qhist.keys():
            if mapq >= histv:
                qhist[histv]= qhist[histv] + mapq_dict[mapq]
    qkeys= sorted([x for x in qhist.keys()])
    lhist= []
    for q in qkeys:
        lhist.append([q, qhist[q]])
    return(lhist)

def nm_count(nm_dict, nm_list):
    """
    nm_dict= {0:10, 1:5, 2:8, 3:4, 5:10, 6:4}
    nm_list= [0, 1, 2, 3, 4]
    nm_count(nm_dict, nm_list)
    >>> [[0, 10], [1, 5], [2, 8], [3, 4], [4, 14]]
    """
    nm_counter= {}
    for k in nm_list:
        nm_counter[k]= 0
    for nm in nm_dict.keys():
        if nm in nm_counter.keys():
            nm_counter[nm]= nm_dict[nm]
        else:
            nm_counter[nm_list[-1]]= nm_counter[nm_list[-1]] + nm_dict[nm]
    nm_counted= []
    for k in nm_list:
        nm_counted.append([k, nm_counter[k]])
    return(nm_counted)    

def hist(mapq_dict, bins= [(0,0), (1,1), (2,2), (3,3), (4,10000)]):
    """
    mapq_dict= {0: 10, 1: 5, 2: 8, 3: 4, 5: 10, 6: 4}
    >>> [[0, 10], [1, 5], [2, 8], [3, 4], [4, 14]]
    """
    xbins= {}
    for k in bins:
        xbins[k]= 0
    for k in mapq_dict.keys():
        for bin in xbins.keys():
            if k >= bin[0] and k <= bin[1]:
                xbins[bin]= xbins[bin] + mapq_dict[k]
    qkeys= sorted([x for x in xbins.keys()])
    lhist= []
    for q in qkeys:
        lhist.append([q[1], xbins[q]])
    lhist[-1][0]= lhist[-2][0] + 1
    return(lhist)

def flatten(lst):
    """
    Returns a flat list from a list containing nested lists
    """
    for elem in lst:
        if type(elem) in (tuple, list):
            for i in flatten(elem):
                yield(i)
        else:
            yield(elem)

def write_stats(fstats, header):
    """
    Read the attributes in object fstats (of class File_Stats) and return them
    in a *list* in the order given in header.
    """
    data_line= []
    keys= header
    for k in keys:
        if type(fstats.__dict__[k]) in (str, int, float) :
            data_line.append(fstats.__dict__[k])
        else:
            data_line= data_line + fstats.__dict__[k] 
    return(data_line)

def lists2rdf(lists):
    """
    Convert a list of lists to a string that R can read into a dataframe.
    The first list is the header.
    """
    bamqc= "bamqc<- data.frame(stringsAsFactors= FALSE, "
    for col in range(0, len(lists[0])):
        for row in range(0, len(lists)):
            if row == 0:
                bamqc= bamqc + str(lists[row][col]) + '= c('
            else:
                if type(lists[row][col]) in (int, float):
                    " Don't quote if element is numeric"
                    bamqc= bamqc + str(lists[row][col]) + ', '
                else:
                    " Quote evrything else "
                    bamqc= bamqc + '"' + str(lists[row][col]) + '", '
        bamqc= bamqc[0:-2]
        bamqc= bamqc + '), '
    bamqc= bamqc[0:-2] + ')'
    return(bamqc)

def meanstdv(x):
    """
    Return mean and sd of a list.
    See http://www.physics.rutgers.edu/~masud/computing/WPark_recipes_in_python.html
    """
    from math import sqrt
    n, mean, std = len(x), 0, 0
    for a in x:
        mean = mean + a
        mean = mean / float(n)
        for a in x:
            std = std + (a - mean)**2
            std = sqrt(std / float(n-1))
    return(mean, std) 

# -----------------------------------------------------------------------------
# Code for quantile found at http://adorio-research.org/wordpress/?p=125
# -----------------------------------------------------------------------------

from math import modf, floor
"""
File	quantile.py
Desc	computes sample quantiles
Author  Ernesto P. Adorio, PhD.
		UPDEPP (U.P. at Clarkfield)
Version 0.0.1 August 7. 2009
"""

def quantile(x, q,  qtype = 7, issorted = False):
	"""
	Args:
	   x - input data
	   q - quantile
	   qtype - algorithm
	   issorted- True if x already sorted.
	Compute quantiles from input array x given q.For median,
	specify q=0.5.
	References:
	   http://reference.wolfram.com/mathematica/ref/Quantile.html
	   http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile
	Author:
	Ernesto P.Adorio Ph.D.
	UP Extension Program in Pampanga, Clark Field.
	"""
	if not issorted:
		y = sorted(x)
	else:
		y = x
	if not (1 <= qtype <= 9):
	   return None  # error!
	# Parameters for the Hyndman and Fan algorithm
	abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
			(0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
			(0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3
			(0,   0, 0, 1), # California linear interpolation, R type 4
			(0.5, 0, 0, 1), # hydrologists method, R type 5
			(0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6
			(1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
			(1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
			(3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
		]
	a, b, c, d = abcd[qtype-1]
	n = len(x)
	g, j = modf( a + (n+b) * q -1)
	if j < 0:
		return y[0]
	elif j >= n:
		return y[n-1]   # oct. 8, 2010 y[n]???!! uncaught  off by 1 error!!!
	j = int(floor(j))
	if g ==  0:
	   return y[j]
	else:
	   return y[j] + (y[j+1]- y[j])* (c + d * g)
def Test():
	x = [11.4, 17.3, 21.3, 25.9, 40.1, 50.5, 60.0, 70.0, 75]
	for qtype in range(1,10):
		print qtype, quantile(x, 0.35, qtype)
#if __name__ == "__main__":
#	Test()
# -----------------------------------------------------------------------------
#                       End code for quantile
# -----------------------------------------------------------------------------
#                       Start processing files
# -----------------------------------------------------------------------------

nm_bins= []
for n in LIMITS_NM[0:-1]:
    nm_bins.append((n, n))
nm_bins.append((LIMITS_NM[-1], 10000))

bamqc_stats= [colnames] ## List (of lists) to hold output

#args.output.write('\t'.join(colnames)+ '\n')

for bam in args.input:
    if outfilename != '<stdout>':
        print('Processing file: %s' %(bam))
    fstats= File_Stats()
    mapq_sample= [] ## Accumulator for sample of mapq scores.
    length_sample= []
    fstats.filename= os.path.split(bam)[1]
    bamfile = pysam.Samfile( bam, "rb" )
    n= 0
    for AlignedRead in bamfile:
        fstats.nreads_tot += 1
        if AlignedRead.flag == 0 or 4 % AlignedRead.flag != 0:
            fstats.nreads_aln += 1
#        else:
#            fstats.nreads_unaln += 1
#            continue
        if AlignedRead.mapq == 255:
            fstats.nreads_mapq_255 = fstats.nreads_mapq_255 + 1 
        else:
            fstats.mapq[AlignedRead.mapq]= fstats.mapq.get(AlignedRead.mapq, 0) + 1
        if n % 100 == 0:
            if AlignedRead.rlen != 0:
                length_sample.append(AlignedRead.rlen)
            if AlignedRead.mapq != 255:
                mapq_sample.append(AlignedRead.mapq)
        try:
            fstats.nreads_nm[AlignedRead.opt('NM')]= fstats.nreads_nm.get(AlignedRead.opt('NM'), 0) + 1
        except:
            fstats.nreads_nm_na= fstats.nreads_nm_na + 1
        n += 1
#        if n >= 10000:
#            break
    fstats.mapq= [ x[1] for x in cumdist(fstats.mapq, LIMITS_MAPQ) ]
    fstats.nreads_nm= [ x[1] for x in hist(fstats.nreads_nm, nm_bins) ]
    fstats.perc_aln= round((float(fstats.nreads_aln) / fstats.nreads_tot )* 100, 2)
    fstats.mapq_quantiles= [ quantile(mapq_sample, q) for q in QUANTILES ]
    fstats.median_length= round(quantile(length_sample, 0.5), 2)
    fstats.len_sd= round(meanstdv(length_sample)[1], 2)
    bamfile.close()
    bamqc_stats.append(write_stats(fstats, header))

for line in bamqc_stats:
    args.output.write('\t'.join([str(x) for x in line]) + '\n')
args.output.close()

if args.nograph is True:
    sys.exit()

# ------------------------------------------------------------------------------
#                            R module for plotting
#
# Input: list of lists with the stats output (bamqc_stats)
# ------------------------------------------------------------------------------

if os.path.split(outfilename)[0] != '':
    os.chdir(os.path.split(outfilename)[0])
if outfilename == '<stdout>':
    outfilename= 'bamqc'
rgraph= """

# ------------------------------------------------------------------------------
# Script for plotting bamqc.py output
# USAGE:
# R CMD BATCH thisfile.R
# ------------------------------------------------------------------------------

#if (is.element("RColorBrewer", installed.packages()) == FALSE){
#     install.packages("RColorBrewer")
#}
library("RColorBrewer")

## Uncomment to read from bamqc.py output
## bamqc<- read.table(file= '...', sep= '\t', header= TRUE)

## Uncomment not to read string from bamqc.py
%(bamqc)s


bamqc$filename<- sub('\\\.bam$', '', bamqc$filename, perl= TRUE) ## Need three backslashes
print(bamqc)
indx<- grep('^mapq_quant\\\.', names(bamqc), perl= TRUE)
bx<- t(bamqc[,indx])
indx.nm<- grep('^nm\\\.', names(bamqc), perl= TRUE)
nm<- t(bamqc[,rev(indx.nm)]) / 1000
print(nm)

WIDTH= 480
HEIGHT= (WIDTH/15) * nrow(bamqc)
MTEXTCEX= 0.8

png('%(pngfile)s', width= WIDTH, height= HEIGHT)

par(mfrow= c(1,3), oma= c(2, 10, 5, 1), mar= c(1,0.5,1,0.5), mgp= c(2.5, 0.5, 0), las= 1)

barplot(t(bamqc[,c("n", "aln", "mapq.15", "mapq.20", "mapq.30")])/1000, horiz= TRUE,
    names.arg= bamqc$filename, beside= TRUE, space=c(-1, 0.5), col= brewer.pal(5, "OrRd"))
axis(side= 3)
mtext(text= 'No. of reads (x1000)\n(all, aligned, mapq 15, 20, 30)', side= 3, line= 2, cex= MTEXTCEX)
grid(ny= NA)
## legend('topleft', legend= c('Tot.', 'Aln.', 'Aln. map <= 15', '<= 20', '<= 30'), horiz= FALSE, inset= c(0, -0.03), bty= 'n', col= col, xpd= NA, pch= 19)
boxplot(bx, yaxt= 'n', horizontal= TRUE, boxvex= 0.2, names= FALSE, ylim= c(0, ifelse(max(bx) > 40, max(bx), 40)))
grid(ny= NA)
axis(side= 3)
mtext(text= 'Quantiles of mapq scores\n(.05, .25, 0.5, .75, .95)', side= 3, line= 2, cex= MTEXTCEX)

barplot(nm, horiz= TRUE, beside= TRUE, space=c(-0.5, 0.5), col= brewer.pal(7, "OrRd"))
grid(ny= NA)
axis(side= 3)
mtext(text= 'No. reads (x1000)\nedit distance: 0, 1 to 5, 6+', side= 3, line= 2, cex= MTEXTCEX)
dev.off()

""" %{'bamqc': lists2rdf(bamqc_stats), 'pngfile': outfilename + '.png'}

tmpfile= outfilename + '.R'
tmp_fh= open(tmpfile, 'w')
tmp_fh.write(rgraph)
tmp_fh.close()

cmd= 'R CMD BATCH %s' %(tmpfile)
p = subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE)
sys.exit()

# ------------------------------------------------------------------------------
#                                   Tritume
# ------------------------------------------------------------------------------

def dict2vect(mydict):
    """ Convert python dict to str formatted as R named vector """
    myk= mydict.keys()
    myk.sort()
    rvect= 'c('
    for k in myk:
        rvect= rvect + '"' + str(k) + '"=' + str(mydict[k]) + ', '
    rvect= rvect.strip(', ') + ')'
    return(rvect)

print('\t'.join(sorted(File_Stats().__dict__)))


