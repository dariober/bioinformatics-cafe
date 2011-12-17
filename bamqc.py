#!/usr/bin/python   

import pysam
import sys

LIMITS_MAPQ= [0, 5, 10, 15, 20, 25, 30, 35, 40]
LIMITS_NM= [0, 1, 2, 3, 4, 5, 6]
bamfiles= sys.argv[1:]

class File_Stats:
    def __init__(self):
        """ Attributes here will be columns in output table """
        self.filename= ''
        self.nreads_tot= 0
        self.nreads_aln= 0
        self.nreads_mapq= 0
        self.nreads_nm= {}
        self.mapq= {}
header=  ['filename', 'nreads_tot', 'nreads_aln', 'mapq', 'nreads_nm']
colnames= ['filename', 'nreads_tot', 'nreads_aln', '\t'.join(['mapq.'+str(x) for x in LIMITS_MAPQ]), '\t'.join(['nm.'+str(x) for x in LIMITS_NM])]

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

def hist(mapq_dict, bins= {(0,0):0, (1,1):0, (2,2):0, (3,3):0, (4,10000):0}):
    """
    mapq_dict= {0: 10, 1: 5, 2: 8, 3: 4, 5: 10, 6: 4}
    >>> [[0, 10], [1, 5], [2, 8], [3, 4], [4, 14]]
    """
    qhist= bins
    for k in mapq_dict.keys():
        for bin in qhist.keys():
            if k >= bin[0] and k <= bin[1]:
                qhist[bin]= qhist[bin] + mapq_dict[k]
    qkeys= sorted([x for x in qhist.keys()])
    lhist= []
    for q in qkeys:
        lhist.append([q[1], qhist[q]])
    lhist[-1][0]= lhist[-2][0] + 1
    return(lhist)

def write_stats(fstats, header):
    data_line= []
    ##keys= sorted(fstats.__dict__.keys())
    keys= header
    for k in keys:
        data_line.append(str(fstats.__dict__[k]))
    print('\t'.join(data_line))

nm_bins= {}
for n in LIMITS_NM[0:-1]:
    nm_bins[(n, n)]= 0
nm_bins[LIMITS_NM[-1], 10000]= 0

print('\t'.join(colnames))
for bamfile in bamfiles:
    fstats= File_Stats()
    fstats.filename= bamfile
    bamfile = pysam.Samfile( bamfile, "rb" )
    n= 0
    for AlignedRead in bamfile:
        fstats.nreads_tot += 1
        if AlignedRead.flag == 0 or 4 % AlignedRead.flag != 0:
            fstats.nreads_aln += 1
#        else:
#            fstats.nreads_unaln += 1
#            continue
        fstats.mapq[AlignedRead.mapq]= fstats.mapq.get(AlignedRead.mapq, 0) + 1
        fstats.nreads_nm[AlignedRead.opt('NM')]= fstats.nreads_nm.get(AlignedRead.opt('NM'), 0) + 1
        n += 1
        if n >= 1000:
            break
    fstats.mapq= '\t'.join( [ str(x[1]) for x in cumdist(fstats.mapq, LIMITS_MAPQ) ] )
    fstats.nreads_nm= '\t'.join( [ str(x[1]) for x in hist(fstats.nreads_nm, nm_bins) ] )
    bamfile.close()
    write_stats(fstats, header)
sys.exit()


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


