import sys
import pybedtools
import os
import tempfile
import subprocess
import pycoverage
import gzip
import re
import pympileup

##COUNT_HEADER= ['depth', 'A', 'a', 'C', 'c', 'G', 'g', 'T', 't', 'N', 'n', 'Z', 'z']
COUNT_HEADER= ['A', 'a', 'C', 'c', 'G', 'g', 'T', 't', 'N', 'n', 'Z', 'z']

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

def mpileup_cmd(bamlist, region, fasta= None, mpileup= 'samtools mpileup'):
    """DEPRECATED: Use mpileup_java_cmd instead()
    Compile a command string to execute samtools mpileup
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
    r= '-r ' + region.chrom + ':' + str(region.start + 1) + '-' + str(region.end)
    if fasta:
        f= '-f %s' %(fasta)
    else:
        f= ''        
    cmd= '%(mpileup)s %(f)s -BQ0 -d10000000 %(r)s %(bamlist)s' %{'mpileup': mpileup, 'f': f, 'r': r, 'bamlist': ' '.join(bamlist)}
    return(cmd)

def mpileup_java_cmd(bamlist, region, fasta= None, mpileup= 'samtools mpileup'):
    """Compile a command string to execute samtools mpileup piped to
    java mpileupParser
    bamlist:
        List of input bams
    region:
        Pybedtools interval to get chrom start, end position
    f:
        FASTA file to get sequence from.
    mpileup:
        String with the full path to mpileup or just samtools mpileup
        if the program is on path
    Return:
        String to pass to subprocess.Popen(). Subprocess will return a string
        formatted as a dict like:
    {'chrom': 'chr7', 'pos': 5566778, 'base': 'N',
    0: {'A': 0, 'a': 0, 'C': 1, 'c': 0, 'G': 0, 'g': 0, 'T': 0, 't': 0, 'N': 0, 'n': 0, 'Z': 1, 'z': 0},
    1: {'A': 0, 'a': 0, 'C': 1, 'c': 0, 'G': 0, 'g': 0, 'T': 0, 't': 0, 'N': 0, 'n': 0, 'Z': 1, 'z': 0},}
    The numeric keys are one for each bamfile passed to mpileup
    """
    r= '-r ' + region.chrom + ':' + str(region.start + 1) + '-' + str(region.end)
    if fasta:
        f= '-f %s' %(fasta)
    else:
        f= ''
    mpileupParserPath= os.path.join(os.path.split(pympileup.__file__)[0], 'java_code/mpileupToNucCounts')
    
    cmd= '%(mpileup)s %(f)s -BQ0 -d10000000 %(r)s %(bamlist)s | java -classpath %(mpileupParserPath)s mpileupParser' %{'mpileup': mpileup,
            'f': f, 'r': r, 'bamlist': ' '.join(bamlist), 'mpileupParserPath': mpileupParserPath}
    return(cmd)

                
def pileupBaseCallsToNucs(bases, refbase):
        """Parses the string of read bases from mpileup output to return the count
        of A, a, C, c, T, t, G, g, N, n and sum of them Z, z. stranded.
        refbase:
            The reference base for this position. This is the 3rd column in mpileup.
        Return:
            Dictionary with where values are counts and keys A, C, G,
            T, N, Z (this latter being the sum).
        
        See also:
            http://samtools.sourceforge.net/pileup.shtml
            Quoting:
            *dot* stands for a match to the reference base on the forward
            strand
            *comma* stands for a match on the reverse strand,
            `ACGTN' for a mismatch on the forward strand and
            `acgtn' for a mismatch on the reverse strand. 
            
            Also at the read base column, a symbol `^' marks the start of a read
            segment which is a contiguous subsequence on the read separated by
            `N/S/H' CIGAR operations. The ASCII of the character following `^'
            minus 33 gives the mapping quality. A symbol `$' marks the end of a
            read segment.
        """
        nuc_counts= {}
        ## Remove the char after '^' since this is mapping quality not base.
        refbase= refbase.upper()
        callDict= {'A': 0, 'a': 0, 'C': 0, 'c': 0, 'G': 0, 'g': 0, 'T': 0, 't': 0, 'N': 0, 'n': 0, '.': 0, ',': 0}
        keys= tuple(callDict)
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
        callDict[refbase] += callDict['.']
        callDict[refbase.lower()] += callDict[',']
        callDict['Z']= sum((callDict[x] for x in ('A', 'C', 'G', 'T', 'N')))
        callDict['z']= sum((callDict[x] for x in ('a', 'c', 'g', 't', 'n')))
        return(callDict)

def parse_pileup(pileup_line, bams):
    """DEPRECATED: Dict is returned by java code (actually string fomratted as dict)
    Parse a pileup line (str) typically returned by pysam.mpileup() or
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

def pileupToBed(pdict, bams, count_header= COUNT_HEADER):
    """Convert the dictionary produced by parse_pileup to a list suitable to
    be written as bedfile. The bed line as:
    <chrom> <pos-1> <pos> <nuc.1> <nuc.2> ... <nuc.n>
    bams:
        list of bam files *in the same order* as passed to samtools.
        Bamfiles are associated to nuc counts using their index in the list. 
    """
    bedlist= [pdict['chrom'], pdict['pos']-1, pdict['pos']]
    bam_idx= range(0, len(bams))
    for c in COUNT_HEADER:
        for idx in bam_idx:
            bedlist.append(pdict[idx][c])
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

def make_dummy_mpileup(chrom, start, end, nbams, count_header= COUNT_HEADER):
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
    zeros= [0] * len(count_header) * nbams
    bedline= [chrom, start, end] + zeros
    return(bedline)

def bamlist_to_mpileup(mpileup_name, mpileup_grp_name, bamlist, region, nwinds, fasta, RPM, regionWindows, samtools, groupFun= 'mean', count_header= COUNT_HEADER):
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
    for h in count_header:
        header.extend([x + '.' + h for x in bamlist])
    header= '\t'.join(header)
    mpileup_bed= open(mpileup_name, 'w')
    if RPM:
        libsizes= getLibrarySizes(bamlist, samtools_path= samtools)
        libsizes= [libsizes[x] for x in bamlist]
    cmd= mpileup_java_cmd(bamlist= bamlist, region= region, fasta= fasta, mpileup= os.path.join(samtools, 'samtools mpileup'))
    proc= subprocess.Popen(cmd, shell= True, stdout=subprocess.PIPE, stderr= subprocess.PIPE)
    nlines= 0
    while True:
        ## Use this while loop to avoid reading in memory all the output of mpileup.
        ## See also http://stackoverflow.com/questions/2804543/read-subprocess-stdout-line-by-line
        line= proc.stdout.readline()
        sys.stdout.flush()
        if not line:
            break
        pd= eval(line)
#        pd= parse_pileup(line, bamlist)
        bedline= pileupToBed(pd, bamlist)
        cnt_indx= 3 ## Column index in bedline where counts start (4th columns)
        if RPM:
            bedline= bedline[0:cnt_indx] + rpm(bedline[cnt_indx:], libsizes)
        mpileup_bed.write('\t'.join([str(x) for x in bedline]) + '\n')
        nlines += 1
    stdout, stderr= proc.communicate()
    if proc.returncode != 0:
        print('\n' + stderr)
        print('samtools exit code: ' + str(proc.returncode) + '\n')
        raise Exception('Failed to execute:\n%s' %(cmd))
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
