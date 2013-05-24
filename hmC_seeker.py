#!/usr/bin/env python

import argparse
import sys
import os
import pybedtools
## from scipy import stats
import fisher
import math

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Compare oxBS and BS bedgraphs to extract hydroxymethylation.
    Input files must be in bedgraph format where the 4th column is the percentage
    methylation
    
    - Complement the two bedgraphs so that all the positions covered in one file are present
      in the other. This way you get the same number and position of windows. 
        - merge bs, oxbs
        - left-outer-join ["merged" left join "BS"] and ["merged" left join "oxBS"]
        
    - Divide bedgraph files in windows and summarize % methylation in windows. Use
      mean or median or sum of reads? (see also compressBedByWindows.py)
    
    - Calculate difference `BS - oxBS`
    
    - Plot distribution of differences 
    
    - Compute distribution of negative percentages (distribution of false positives).
      Get the 95% quantile x95. Everything larger than abs(x95) will be considered
      true positive hmC    

INPUT FORMAT
    Bedgraph with additional columns as produced by mpileup2methylation.py additional
    columns (e.g. strand are ignored):
    
    <chrom>  <start>  <end>  <pct met (unusued)>  <cnt methylated>  <tot count>

OUTPUT
    Bed file with columns
    
    <chrom> <w/start> <w/end> <BS cnt met> <BS cnt tot> <BS pct met> <oxBS cnt met> <oxBS cnt tot> <oxBS pct met> <fisher left tail p> <<fisher right tail p>>
    
    "fisher left tail p": pvalue for H "methylation % BS <= % oxBS"
    "fisher right tail p": pvalue for H "methylation % BS >= % oxBS" (p~0 means there is hmC)
    
    Example:
    chrM    1       11      0       31      0.0     0       36      0.0     1.0     1.0
    chrM    11      21      0       111     0.0     1       113     0.8772  0.5067  1.0
    chrM    21      31      0       206     0.0     0       230     0.0     1.0     1.0
    chrM    31      41      0       255     0.0     1       348     0.2865  0.5778  1.0

    
EXAMPLE    

    
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--bs',
                   required= True,
                   help='''Input bedgraph of the BS library

''')

parser.add_argument('--oxbs',
                   required= True,
                   help='''Input bedgraph of the oxBS library

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

parser.add_argument('--fisher', '-f',
                    action= 'store_true',
                   help='''Execute Fisher exact test.

''')

args= parser.parse_args()

def reduceColumns(feature):
    """Remove redundant columns from BedTool feature and replace missing values '.'
    with '0'.
    """
    newfeature= pybedtools.Interval(feature.chrom, feature.start, feature.end, name= 'NA')
    if feature[3] == '.':
        newfeature.append('\t'.join(['0', '0', '0']))
    else:
        newfeature.append('\t'.join(feature[6:9]))    ## This is where input columns are added
    return(newfeature)
   
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
    summWinds= intsct.groupby(g= [1,2,3], c= [10, 11, 12], o= [ops, ops, ops], stream= False)
    return(summWinds)

def pctMet(M, tot):
    """Return % methylated `M / (M+m)` or NA if 0/0.
    M, m: int for count methylated and unmethylated
    """
    if tot == 0:
        return('NA')
    elif M > tot:
        sys.exit('Count methylated > Total coverage: %s, %s' %(M, tot))
    else:
        pct= round(100*(float(M) / tot), 4)
        return(pct)
    
    
def fisherExact(line, idx):
    """Apply fisher exact test to appropriate columns of bed line. Columns are
    selected with the indexes in idx
    """
    ## cnt= [[int(line[3]), int(line[4])], [int(line[5]), int(line[6])]]
    fet= fisher.pvalue(int(line[idx[0]]), int(line[idx[1]]), int(line[idx[2]]), int(line[idx[3]]))
    pvalues= [str(round(fet.left_tail, 4)), str(round(fet.right_tail, 4))]
    line.append('\t'.join(pvalues))
    return(line)
    
# -----------------------------------------------------------------------------

## Column indexes to use
INFMT= {'chrom': 0, 'start': 1, 'end': 2, 'score': 3, 'cnt_M': 4, 'cnt_tot': 5}

colIdx= sorted(INFMT.values())

bs_bdg= pybedtools.BedTool(args.bs)
bs_bdg= pybedtools.BedTool( [feature[i] for i in colIdx] for feature in bs_bdg).saveas()
oxbs_bdg= pybedtools.BedTool(args.oxbs)
oxbs_bdg= pybedtools.BedTool( [feature[i] for i in colIdx] for feature in oxbs_bdg).saveas()

## Get the union of all positions.
mergeBS= bs_bdg.cat(oxbs_bdg, postmerge=True, d= -1)

## Add to each lib the missing positions
bs_bdg_loj= mergeBS.intersect(bs_bdg, wb= True, loj= True)
bs_bdg_loj= bs_bdg_loj.each(reduceColumns).saveas()

oxbs_bdg_loj= mergeBS.intersect(oxbs_bdg, wb= True, loj= True)
oxbs_bdg_loj= oxbs_bdg_loj.each(reduceColumns).saveas()

BSwindows= compressBed(bs_bdg_loj, window_size= args.window_size, step_size= args.step_size, ops= args.ops)
oxBSwindows= compressBed(oxbs_bdg_loj, window_size= args.window_size, step_size= args.step_size, ops= args.ops)

windowsPaste= BSwindows.intersect(b= oxBSwindows, wa= True, wb= True)

#SY005_2x2_hmc_Q38_indexed       0       1       1.25229999999999996874  323     25793   SY005_2x2_hmc_Q38_indexed       0       1       1.6546000000000000707   142     8582
#SY005_2x2_hmc_Q38_indexed       2       3       98.7604999999999932925  25974   26300   SY005_2x2_hmc_Q38_indexed       2       3       98.8041999999999944748  8676    8781
#SY005_2x2_hmc_Q38_indexed       4       5       1.11739999999999994884  294     26312   SY005_2x2_hmc_Q38_indexed       4       5       1.18359999999999998543  104     8787
#SY005_2x2_hmc_Q38_indexed       5       6       97.5794999999999959073  25680   26317   SY005_2x2_hmc_Q38_indexed       5       6       97.0776000000000038881  8537    8794

windowsPaste= pybedtools.BedTool((f[INFMT['chrom']],
                                  f[INFMT['start']],
                                  f[INFMT['end']],
                                  f[INFMT['cnt_M']],
                                  int(f[INFMT['cnt_tot']]) - int(f[INFMT['cnt_M']]),   ## Count unmet BS
                                  pctMet(int(f[INFMT['cnt_M']]), int(f[INFMT['cnt_tot']])),   ## Pct methylated BS
                                  f[len(colIdx) + INFMT['cnt_M']],
                                  int(f[len(colIdx) + INFMT['cnt_tot']]) - int(f[len(colIdx) + INFMT['cnt_M']]), ## Count unmet oxBS
                                  pctMet(int(f[len(colIdx) + INFMT['cnt_M']]), int(f[len(colIdx) + INFMT['cnt_tot']])) ## Pct methylated oxBS                                
                                  ) 
                                  for f in windowsPaste)

if args.fisher:
    windowsPaste= windowsPaste.each(fisherExact, [3,4, 6,7]).saveas()
for line in windowsPaste:
    sys.stdout.write(str(line))

sys.exit()
