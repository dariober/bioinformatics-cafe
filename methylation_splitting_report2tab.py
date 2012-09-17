#!/home/berald01/.local/bin/python

import argparse
import sys
import os
import re
import inspect
import traceback

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Convert the report from "methylation_report --report" to tabular format, tab separated, suitable
    for import to spreadsheet or database.
    
EXAMPLE
    ## Parse all repports ending in clean_splitting_report.txt, output header line
    ls *clean_splitting_report.txt | methylation_splitting_report2tab.py --first_header -
TODO:
    - Add parsing to identify strand-specificity, single/paired end
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('infile',
                   nargs= '+',
                   help='''One or more reports to parse. Use - to read the list of files from stdin
                   ''')

parser.add_argument('--columns',
                   action= 'store_true',
                   help='''Print report as two columns: Headers, Values
                   ''')

parser.add_argument('--first_header',
                   action= 'store_true',
                   help='''The first line to be printed out is the header.
                   ''')

args = parser.parse_args()

# ------------------------------------------------------------------------------

# The functions get_... read the report and return the value for the table header
# following get_ (e.g. get_bismark_report_for() returns the value for bismark_report_for)
# ARG report_list is the report as list as produced by (no leading and trailing blanks, no empty lines):
# report_list= open(bismark_report).readlines()
# report_list= [x.strip() for x in report_list if not x.strip() == '']
# Return a tuple as: (header: value) e.g. ('version', 'v0.7.6') 
def get_splitting_report_for(report_list):
    line= report_list[0]
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_no_call_strings(report_list):
    tag= 'Total number of methylation call strings processed: '
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '')
    line= re.sub(' .*', '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_tot_cs(report_list):
    tag= "Total number of C's analysed:\t"
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_mC_cpg(report_list):
    tag= "Total methylated C's in CpG context:\t"
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_mC_chg(report_list):
    tag= "Total methylated C's in CHG context:\t"
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_mC_chh(report_list):
    tag= "Total methylated C's in CHH context:\t"
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_c2t_cpg(report_list):
    tag= "Total C to T conversions in CpG context:\t"
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_c2t_chg(report_list):
    tag= "Total C to T conversions in CHG context:\t"
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_c2t_chh(report_list):
    tag= "Total C to T conversions in CHH context:\t"
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_perc_mC_cpg(report_list):
    tag= "C methylated in CpG context:\t"
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        return(re.sub('^get_', '',  inspect.stack()[0][3]), 'NA')
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line).rstrip('%')
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_perc_mC_chg(report_list):
    tag= "C methylated in CHG context:\t"
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        return(re.sub('^get_', '',  inspect.stack()[0][3]), 'NA')
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line).rstrip('%')
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_perc_mC_chh(report_list):
    tag= "C methylated in CHH context:\t"
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        return(re.sub('^get_', '',  inspect.stack()[0][3]), 'NA')
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line).rstrip('%')
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

# ------------------------------------------------------------------------------

if args.infile == ['-']:
    args.infile= sys.stdin.readlines()
    args.infile= [x.strip() for x in args.infile]

data= [] ## Store each report a string here as list
for bismark_report in args.infile:
    report_list= open(bismark_report).readlines()
    report_list= [x.strip() for x in report_list if not x.strip() == '']
    try:
        values= [## Each item is tuple (header_name, value)
                get_splitting_report_for(report_list),
                get_no_call_strings(report_list),
                get_tot_cs(report_list),
                get_mC_cpg(report_list),
                get_mC_chg(report_list),
                get_mC_chh(report_list),
                get_c2t_cpg(report_list),
                get_c2t_chg(report_list),
                get_c2t_chh(report_list),
                get_perc_mC_cpg(report_list),
                get_perc_mC_chg(report_list),
                get_perc_mC_chh(report_list),
                ('report_name', bismark_report)]
    except:
        stack_trace = traceback.format_exc()
        print('\nReport "%s" raised an exception:\n' %(bismark_report))
        sys.exit(stack_trace)
        
    data.append(values)
if args.first_header:
    header= [x[0] for x in data[0]]
    print('\t'.join(header))

for report in data:
    if args.columns:
        print('-'*60)
        for v in report:
            print(v[0] + ':\t' + str(v[1]))
    else:
        print('\t'.join([str(v[1]) for v in report]))

sys.exit()

