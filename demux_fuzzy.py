#!/home/berald01/.local/bin/python

import Levenshtein
import sys
import argparse
import gzip
import subprocess
import operator
try:
    import sblab
except:
    pass
    ## raise Warning('Module sblab not found on this system.')

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    De-multiplex a FASTQ file into separate files on the bases of the barcode sequence
    found on the read name. The output files are gzipped.
    
    The barcode is extracted from the header line of each sequence by extracting the substring
    between rightmost '#' and '/'.
    From this substring, the first 6 bases are extracted as the barcode.
    barcode_sequence= hline[(hline.rfind('#')+1) : hline.rfind('/')][0:6]
    
    demux_fuzzy.py rescues imperfect matches provided that
    1) A unique best match can be found
    2) The edit distance between sequence and barcodes is less than a given threshold
       (1 by default).

    The edit distance is the Levenshtein distance as implemented in the python package
    Levenshtein
    
    The sample sheet gives the link between barcode and output (like demux).
    It has no header and the first two columns (1st barcode, 2nd file name) must be separated
    by space (*not* TAB). Additional columns are ignored.
    
-------sample sheet example ---------
TAGCTTA demultiplex_file-1.fq
ACTTGAA demultiplex_file-2.fq
-------------------------------------
    
USAGE:
    demux_fuzzy.py -f <fastqfile> -s <sample sheet>
    
    ## Read from stdin a zipped file. Use -f -:
    gunzip -c fastq.fq.gz | demux_fuzzy.py -f - -s samplesheet.txt

    ## Upload report to sblab. Report will be in file mysheet.demux_report
    demux_fuzzy.py -f myfastq.fq.gz -s mysheet.txt --pgupload --report
    
    ## Only upload a previoiusly produced report. Options -f and -s are ignored but required (bug).
    demux_fuzzy.py -f bla -s bla -u test.demux.txt.demux_report

REQUIREMENTS:
    Python with package Levenshtein (http://pypi.python.org/pypi/python-Levenshtein/)
    and argparse (http://pypi.python.org/pypi/argparse)

TODO
   - Allow to output in uncompressed format.
   - Allow for insertions/deletions by searching the whole barcode sequence in the header?
     (Probably not worth the effort)
     ...
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--fastq', '-f',
                   type= str,
                   help='''Input FASTQ file to de-multiplex. Use - to read from stdin
                                   
                   ''',
                   required= True)

parser.add_argument('--samplesheet', '-s',
                   type= str,
                   help='''Sample sheet with barcodes (1st colulmn) and output
files (2nd column). Additional columns are ignored.
Memo: Columns are space (not tab) separated.
                                   
                   ''',
                   required= True)

parser.add_argument('--distance', '-d',
                   type= int,
                   help='''Maximum edit distance between barcodes found in the sample
sheet and sequence read from the fastq file. Default 1. (Memo: ambiguos barcodes
are always discarded)
                                   
                   ''',
                   required= False,
                   default= 1)

parser.add_argument('--report', '-r',
                   type= str,
                   nargs= '?',
                   help='''Report file where to write some QC.
Use - to send to stdout (default). Leave blank to send to file <sample sheet>.demux_report
                   ''',
                   required= False,
                   default= '-')

parser.add_argument('--pgupload', '-p',
                   action= 'store_true',
                   help='''Upload report to postgres sblab.main.demux_report
                   ''')

parser.add_argument('--pgupload_only', '-u',
                   type= str,
                   required= False,
                   help='''Only upload a report file. Using sblab module.
                   ''')

args = parser.parse_args()

# -----------------------------------------------------------------------------

barcode_dict= { "ACAGTG": "TruSeq-5",
                "ACTGAT": "TruSeq-25",
                "ACTTGA": "TruSeq-8",
                "AGTCAA": "TruSeq-13",
                "AGTTCC": "TruSeq-14",
                "ATCACG": "TruSeq-1",
                "ATGTCA": "TruSeq-15",
                "ATTCCT": "TruSeq-27",
                "CAGATC": "TruSeq-7",
                "CCGTCC": "TruSeq-16",
                "CGATGT": "TruSeq-2",
                "CGTACG": "TruSeq-22",
                "CTTGTA": "TruSeq-12",
                "GAGTGG": "TruSeq-23",
                "GATCAG": "TruSeq-9",
                "GCCAAT": "TruSeq-6",
                "GGCTAC": "TruSeq-11",
                "GTCCGC": "TruSeq-18",
                "GTGAAA": "TruSeq-19",
                "GTGGCC": "TruSeq-20",
                "GTTTCG": "TruSeq-21",
                "TAGCTT": "TruSeq-10",
                "TGACCA": "TruSeq-4",
                "TTAGGC": "TruSeq-3"}

if args.pgupload_only:
    reportname= args.pgupload_only
    sblab.uplod_demux_fuzzy_report(reportname)
    sys.exit()
if args.report is None:
    reportname= args.samplesheet + '.demux_report'
    fhreport= open(reportname,  'w')
elif args.report == '-':
    fhreport= sys.stdout
else:
    reportname= args.report
    fhreport= open(reportname, 'w')

#------------------------[ Functions ]-----------------------------------------

def read_fastq_line(fastq_fh):
    """
    Given the file handle fastq_fh (e.g. fastq_fh= open('myfastqfile')), reads
    a chunk of 4 lines.
    """
    fqline= [fastq_fh.readline().strip()]
    fqline.append(fastq_fh.readline().strip())
    fqline.append(fastq_fh.readline().strip())
    fqline.append(fastq_fh.readline().strip())
    return(fqline)
    
def read_samplesheet(sample_sheet):
    """
    Read the file sample_sheet containing. Format is space or tab separated with barcode
    sequence (first column) and output file (second column). Additional columns
    ignored.
    Returns a dictionary of codes and files:
    {'ACTACT': 'file1.fq', 'ACGACG': 'file2.fq', ...}
    """
    ssheet= open(sample_sheet).readlines()
    ssheet= [x.strip() for x in ssheet if x.strip() != ''] ## Strip leading and trailing blanks and blank lines
    ssheet= [x.split(' ') for x in ssheet]
    codes= [x[0] for x in ssheet]
    if len(codes) != len(set(codes)):
        sys.exit('Duplicate barcodes found in sample sheet')
    code_dict= {}
    print('\nSample sheet: %s' %(sample_sheet))
    for line in ssheet:
        print('\t'.join(line))
        ## For each barcode/file store the following
        if line[1].endswith('.gz'):
            filename= line[1]
        else:
            filename= line[1] + '.gz'
        if len(line[0]) == 7:
            index_seq= line[0][0:6]
        elif len(line[0]) == 6:
            index_seq= line[0]
        else:
            sys.exit('Unexpected barcode: %s' %(line[0]))

        code_dict[index_seq]= [
                  filename,                  ## Output file name
                  gzip.open(filename, 'wb'),  ## Output file handle
                  0                             ## Counter for number of reads in this file
                  ]
    print('')
    return(code_dict)

def spurious_hits(barcode_dict_matches, exclude_codes, totreads= None):
    """
    Returns a lst of top spurious barcodes.
    
    barcode_dict_matches: Dictionary of barcodes and their hit count.
        Like {'ACTGAC': ['TruSeq-X', 10], }
    exclude_codes: List of barcodes expected to be present in the master file
        and therefore don't consider them as spurious.
    totreads: Total reads to compute percent spurious. If None, percent is not
              calculated
    
    OUTPUT: List of list in descending order: [['TruSeq-A', 100], ['TruSeq-B', 90], ... ]
    """
    spurhits= {}
    for k in barcode_dict_matches:
        if k in exclude_codes:
            pass
        else:
            spurhits[k]= barcode_dict_matches[k][::-1] ## Make count to be the first
    sorted_x = sorted(spurhits.iteritems(), key=operator.itemgetter(1), reverse= True)
    sorted_x= [x[1][::-1] for x in sorted_x] ## MAke code name to be first
    if totreads is not None:
        for i in range(0, len(sorted_x)):
            x= sorted_x[i]
            x.append(round(float(x[1]) / float(totreads), 2))
            sorted_x[i]= x
    return(sorted_x)

   
# -----------------------------------------------------------------------------

barcode_dict_matches= {}
for k in barcode_dict:
    barcode_dict_matches[k]= [barcode_dict[k], 0]
barcode_dict_matches['N']= ['Barcodes_with_N', 0]
barcode_dict_matches['no match']= ['No_match', 0]

if args.fastq == '-':
    fh= sys.stdin
elif args.fastq.endswith('.gz'):
    fh= gzip.open(args.fastq)
else:
    fh= open(args.fastq)

code_dict= read_samplesheet(args.samplesheet)
codes= []
for k in code_dict:
    codes.append(k)
    if k[0:6] not in barcode_dict.keys():
        print('WARNING: barcode sequence %s is not in current dictionary' %(k))
n= 0
n_lost= 0
while True:
    fqline= read_fastq_line(fh)
    if fqline[0] == '':
        break
    n += 1
    if n % 1000000 == 0:
        print(n)
    hline= fqline[0]
    
    bcode= hline[(hline.rfind('#')+1) : hline.rfind('/')][0:6]## fqline[0][-9:-3]
    ## Check match between barcode and barcode dict. This is only for reporting/QC
    if bcode in barcode_dict_matches.keys():
        barcode_dict_matches[bcode][1] += 1
    elif 'N' in bcode.upper():
        barcode_dict_matches['N'][1] += 1
    else:
        barcode_dict_matches['no match'][1] += 1
    ## Do the actual demultiplexing
    if bcode in code_dict:
        ## Test for perfect match
        code_dict[bcode][1].write('\n'.join(fqline) + '\n')
        code_dict[bcode][2] += 1
    elif args.distance == 0 or bcode.count('N') > args.distance:
        ## If there are more Ns than allowed by args.distance, don't go to Levenshtein distance at all 
        n_lost += 1
        continue
    else:
        dists= [Levenshtein.distance(bcode, x) for x in codes]
        best_dist= min(dists)
        if best_dist > args.distance or len([x for x in dists if x == best_dist]) > 1:
            n_lost += 1
        else:
            bcode= codes[dists.index(best_dist)]
            code_dict[bcode][1].write('\n'.join(fqline) + '\n')
            code_dict[bcode][2] += 1

fh.close()

## Convert list of lists to string formatted like a postgres array
match_report= sorted(barcode_dict_matches.values())
match_report_array= []
for x in match_report:
    match_report_array.append("{'" + x[0] + "'," + str(x[1]) + '}')
match_report_array='{' + ','.join(match_report_array) + '}'

spur_report= spurious_hits(barcode_dict_matches, codes, n)
spur_report_array= []
for x in spur_report:
    spur_report_array.append("{'" + x[0] + "'," + str(x[1]) + ',' + str(x[2]) + '}')
spur_report_array='{' + ','.join(spur_report_array) + '}'

perc_lost= round(100*(n_lost/float(n)),2)
print('\nTotal reads: %s' %(n))
print('Lost:        %s (%s%%)' %(n_lost, perc_lost))
print('\nReads in:')
for k in code_dict:
    code_dict[k][1].close()
    p= subprocess.Popen('get_file_stats2.py -i %s --md5sum' %(code_dict[k][0]), stdout= subprocess.PIPE, shell= True)
    fstats= p.stdout.read()
    fstats= eval(fstats)
    perc= round(100*(code_dict[k][2]/float(n)),2)
    print('%s\t%s\t%s%%' %(code_dict[k][0], code_dict[k][2], perc,))
    ## Row to print in report 
    reportline= []
    reportline.append(code_dict[k][0]) ## Demultiplexed file name
    reportline.append(args.fastq) ## Master file name
    reportline.append(k) ## Barcode seq
    reportline.append(barcode_dict[k]) ## Barcode name (e.g. TruSeq-1)
    reportline.append(code_dict[k][2]) ## N. reads
    reportline.append(n) ## Tot reads in master file
    reportline.append(perc) ## % in this file
    reportline.append(n_lost) ## Reads lost
    reportline.append(perc_lost) ## Reads lost
    reportline.append(fstats['md5sum'])
    reportline.append(fstats['fsize'])
    reportline.append(fstats['mtime'])
    reportline.append(fstats['path'])
    reportline.append(match_report_array)
    reportline.append(spur_report_array)
    fhreport.write('\t'.join([str(x) for x in reportline]) + '\n')

print('')
fhreport.close()

if args.pgupload and args.report != '-':
    sblab.uplod_demux_fuzzy_report(reportname)
sys.exit()