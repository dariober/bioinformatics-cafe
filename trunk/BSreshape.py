#!/usr/bin/env python

import argparse
import sys
import pybedtools
import re
import os
import tempfile
import glob

DEFAULT= {'cnt_M': 4, 'cnt_tot': 5}

parser = argparse.ArgumentParser(description= """
DESCRIPTION

* Merge all input bedgraph to get "consensus" bedgraph (union)
* Add to each bedgraph the missing positions. Missing positions get 0 counts
* Compress bedgraphs by windows
* Concatenate bedgraphs and reshape to have <chrom> <start> <end> <locus name (expt)> <cnt_M> <cnt_c> <lib type (BS/oxBS)> <library_id>
* Feed cat bedgraph to glmBS()

EXAMPLE
    BSreshape.py -w 10 --oxbs mjb050_E14oxBSAD01.20130301.bs4.bedGraph --bs mjb053_E14BSAD04.20130301.bs4.bedGraph
    
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--bs',
                   required= True,
                   nargs= '+',
                   help='''List of bedgraph files for BS libraries. Wild card
chars will expanded.

''')

parser.add_argument('--oxbs',
                   required= True,
                   nargs= '+',
                   help='''List of bedgraph files for oxBS libraries. Wild card
chars will expanded.

''')

parser.add_argument('--strip', '-S',
                   required= False,
                   help='''Strip this regex from file names to generate library_ids
Default is to use the filename dir. NB: library_id's must be unique.
''')

parser.add_argument('--window_size', '-w',
                   required= True,
                   type= int,
                   help='''Window size to group features. This option passed to
bedtools makewindows.

''')

parser.add_argument('--step_size', '-s',
                   required= False,
                   default= None,
                   type= int,
                   help='''Step size to slide windows. This option passed to
bedtools makewindows. Default step_size= window_size (non-sliding windows)

''')

parser.add_argument('--ops', '-o',
                   required= False,
                   default= 'sum',
                   type= str,
                   help='''Operation to apply to the score column. Default: sum
This option passed to bedtools groupBy.

''')

parser.add_argument('--cnt_M', '-M',
                   required= False,
                   default= DEFAULT['cnt_M'],
                   type= int,
                   help='''Column position for the field "Count methylated".
Default %s. This is 0-based (1st column has index 0).

''' %(DEFAULT['cnt_M']))

parser.add_argument('--cnt_tot', '-t',
                   required= False,
                   default= DEFAULT['cnt_tot'],
                   type= int,
                   help='''Column position for the field "total reads (methylated + unmethylated)".
Default %s. This is 0-based (1st column has index 0).

''' %(DEFAULT['cnt_tot']))

args= parser.parse_args()

if args.cnt_M == args.cnt_tot:
    sys.exit('Columns cnt_M and cnt_tot cannot be the same!')

# -----------------------------------------------------------------------------
class Bedgraph:
    def __init__(self):
        self.filename= None
        self.library_id= None
        self.bsType= None
        self.bdgComplete= None
        self.gffComplete= None
        self.gffWindow= None
        self.bdgGLM= None ## Bed file suitable for glmBS()

def bed2gff(bedFeature, source= None, featureName= None, frame= '.', **kwargs):
    """Convert bed format to gff. This way you are not bound to fixed column order for
    methylation counts.
    Use source BS or oxBS to identify the library type.
    bedFeature:
        BED line from BedTool object
    source:
        Use this string for the source field (2nd column in GFF)
    featureName:
        Use this string for the feature field (3rd column in GFF). Defualt is to
        create the name using <chrom>_<start>_<end>
    frame:
        String for the frame field
    **kwargs:
        Each additional argument becomes an attribute of the GFF line.

    See also
        For the gff format see http://www.sanger.ac.uk/resources/software/gff/spec.html
        pybedtools intervals: http://pythonhosted.org/pybedtools/intervals.html
    """
    attr=[]
    for k in kwargs:
        "Convert kwargs dict to attribute string"
        x= str(k) + '=' + str(kwargs[k])
        attr.append(x)
    attr= '; '.join(attr) + ';'
    if source is None:
        source= '.'
    if featureName is None:
        featureName= '_'.join([bedFeature.chrom, str(bedFeature.start), str(bedFeature.end)])
    if bedFeature.strand == '':
        strand= '.'
    else:
        strand= bedFeature.strand
    gff= [bedFeature.chrom, source, featureName, str(bedFeature.start + 1), str(bedFeature.end), bedFeature.score, strand, frame, attr]
    intvbed= pybedtools.create_interval_from_list(gff)
    return(intvbed)

def BEDfile2GFFfile(bedfile, source= None, featureName= None, frame= '.', attrIdx= None, attrStr= None):
    """Convert BED file object to GFF file object.
    bedfile:
        Input bed file obejct from pybedtools.BedTool()
    featureName:
        If None, the name (3rd col) will be <chrom>_<bed start>_<bed end> otherwise
        it will be fetched from the bed file.
    attrIdx:
        Dictionary in the form {'attr': i}. The key will be an attr of the GFF
        with the value in the column i (0-based) of the bed file
    attrStr:
        Dictionary in the form {'attr': x}. The key will be an attr of the GFF
        with value equal to x
    
    Return a pybedtools.BedTool() file object
    
    E.g.
    bedfile:
    chrX 0 1 100 120
    BEDfile2GFFfile(bedfile, attrIdx= {'cnt_M': 3, 'cnt_tot': 4}, attrStr= {'libType': 'BS', 'libID': 'libX'})
    chrM 
    """
    keyCollision= [k for k in attrIdx if k in attrStr]
    if not keyCollision == []:
        sys.exit('Key collision')
    tmpGff= tempfile.NamedTemporaryFile(prefix= 'tmpGff_BSreshape_', delete= False) ## open('tmp.gff', 'w') 
    for feature in bedfile:
        if source is None:
            source= '.'
        if featureName is None:
            name= '_'.join([feature.chrom, str(feature.start), str(feature.end)])
        else:
            name= feature.name
        if feature.strand == '':
            strand= '.'
        else:
            strand= feature.strand
        attrDict= attrStr
        if feature.score == '':
            score= '.'
        else:
            score= feature.score
        for attr in attrIdx:
            "Collect all the attributes"
            i= attrIdx[attr]
            value= feature[i]
            attrDict[attr]= value
        attr= []
        for k in attrDict:
            "Convert dict dict to attribute string"
            x= str(k) + '=' + str(attrDict[k])
            attr.append(x)
        attr= '; '.join(attr) + ';'
        gff= [feature.chrom, source, name, str(feature.start + 1), str(feature.end), str(score), strand, frame, attr]
        tmpGff.write('\t'.join(gff) + '\n')
    tmpGff.close()
    gffFile= pybedtools.BedTool(tmpGff.name).saveas()
    os.remove(tmpGff.name)
    return(gffFile)
        
def addMissingAttr(feature):
    """Add missing attribute to BD3-GFF line.    
    `feature` is a pybedtools.Interval where first three fields are BED3 format
    and the reamining fileds are GFF.
    Return pybedtools.Interval formatted as GFF
    """
    if feature[4] == '.':
        gffLine= bed2gff(feature, source= bdgObj.bsType, cnt_M= 0, cnt_tot= 0, library_id= bdgObj.library_id, bsType=bdgObj.bsType)
    else:
        gffLine= pybedtools.create_interval_from_list(feature[3:])
    return(gffLine)

def compressBed(inbed, window_size, step_size, ops):
    """    
    """
    ## 1. Get extremes of each chrom
    grp= inbed.groupby(g= [1], c= [2,3], ops= ['min', 'max'], stream= False)
    ## 2. Divide each chrom in windows
    windows= grp.window_maker(b= grp.fn, w= window_size, s= step_size, stream= False)

    ## 3. Assign bed features to windows
    intsct= windows.intersect(b= inbed, wa= True, wb= True, stream= False)
    
    ## 4. Summarize windows
    summWinds= intsct.groupby(g= [1,2,3], c= [7,8], o= [ops, ops], stream= False)
    return(summWinds)

# -----------------------------------------------------------------------------
# Get files and some checks on them
ff= []
for f in args.bs:
    flist= glob.glob(f)
    ff.extend(flist)
bsFiles= sorted(set(ff))
ff= []
for f in args.oxbs:
    flist= glob.glob(f)
    ff.extend(flist)
oxbsFiles= sorted(set(ff))

dupFiles= [f for f in bsFiles if f in oxbsFiles]
if not dupFiles == []: 
    sys.exit('Some files in --oxbs list are also in the --bs list.')
bsfiles= bsFiles+ oxbsFiles
# -----------------------------------------------------------------------------

# Merge all input bedgraphs to get "consensus" bedgraph (union)
unionBed= pybedtools.BedTool(bsfiles[0])
unionBed= unionBed.cat(*bsfiles, postmerge= True)

# Collect info on the libraries
bdgDict={} ## This dict will be {'library_id': Bedgraph()}. It's not really used but keep for future devel
print('\t'.join(['locus', 'library_id', 'bs', 'cnt_M', 'cnt_c']))
for bdg in bsfiles:
    bdgObj= Bedgraph()
    bdgObj.filename= bdg
    bdgObj.library_id= os.path.split(bdgObj.filename)[1]
    if args.strip:
        bdgObj.library_id= re.sub(args.strip, '', bdgObj.library_id)
    if bdgObj.library_id in bdgDict:
        sys.exit('\nFound duplicate file names or library IDs:\n%s\n' %(bdgObj.library_id))
    else:
        bdgDict[bdgObj.library_id]= bdgObj
    if bdg in bsFiles:
        bdgObj.bsType= 'BS'
    else:
        bdgObj.bsType= 'oxBS'
    ## Convert bed to gff
    bed= pybedtools.BedTool(bdgObj.filename)
    bdgObj.gffComplete= BEDfile2GFFfile(bed, source= bdgObj.bsType,
        attrStr= {'library_id': bdgObj.library_id, 'bsType': bdgObj.bsType},
        attrIdx= {'cnt_M': 4, 'cnt_tot': 5})
    ## Add missing positions. After this gffComplete is not a real GFF!
    bdgObj.gffComplete= unionBed.intersect(bdgObj.gffComplete, wb= True, loj= True).saveas()
    ## Assign missing attributes. Now gffComplete is real GFF again
    bdgObj.gffComplete= bdgObj.gffComplete.each(addMissingAttr)
    ## Compress GFF. Need to convert to BED in order to use groupBy
    bed= [[f.chrom, f.start-1, f.end, f.attrs['cnt_M'], f.attrs['cnt_tot']] for f in bdgObj.gffComplete]
    xbed= pybedtools.BedTool(bed).saveas()
    windBed= compressBed(xbed, args.window_size, args.step_size, args.ops)
    ## Convert BED 2 GFF
    bdgObj.gffWindow= BEDfile2GFFfile(windBed, source= bdgObj.bsType,
        attrStr= {'library_id': bdgObj.library_id, 'bsType': bdgObj.bsType},
        attrIdx= {'cnt_M': 3, 'cnt_tot': 4})
    # -------------------------------------------------------------------------
    # Print out file suitable for R/glmBS
    # Column headers must be: library_id, locus, cnt_M, cnt_c, bs
    for feature in bdgObj.gffWindow:
        cnt_c= str(int(feature.attrs['cnt_tot']) - int(feature.attrs['cnt_M']))
        line= [feature[2], feature.attrs['library_id'], feature.attrs['bsType'], feature.attrs['cnt_M'], cnt_c]
        print('\t'.join(line))
sys.exit()

