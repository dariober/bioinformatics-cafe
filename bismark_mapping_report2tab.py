#!/home/berald01/.local/bin/python

import argparse
import sys
import os
import re
import inspect
import traceback

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Convert the Bismark mapping report to tabular format, tab separated, suitable
    for import to spreadsheet or database.
    Currently tested on (some) reports from bismark v0.7.4 and v0.7.6 PE and SE.
    
EXAMPLE
    bismark_mapping_report2tab.py bis_map_rep.txt bis_map_rep2.txt
    
    ## Convert a bunch of reports read by stdin, add header line
    ls *mapping_report.txt | bismark_mapping_report2tab.py --first_header - 
    
    ## For on-screen ease of reading use --columns:
    ls *mapping_report.txt | bismark_mapping_report2tab.py --columns -
    
TODO:

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
def get_bismark_report_for(report_list):
    tag= 'Bismark report for: '
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    line= re.sub('( \(version\: ).*\)$', '', line)
    line= line.split(' and ')
    line= [os.path.split(x)[1] for x in line]
    if len(line) == 1:
        line= line[0]
    elif len(line) == 2:
        line= ', '.join(line)
    else:
        sys.exit('%s: Too many items files found: Should be either 1 or 2 for %s' %(inspect.stack()[0][3], line))
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_full_path(report_list):
    tag= 'Bismark report for: '
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    line= re.sub('( \(version\: ).*\)$', '', line)
    line= line.split(' and ')
    if len(line) == 1:
        line= line[0]
    elif len(line) == 2:
        line= ', '.join(line)
    else:
        sys.exit('%s: Too many items files found: Should be either 1 or 2 for %s' %(inspect.stack()[0][3], line))
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_bismark_version(report_list):
    tag= 'Bismark report for: '
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    line= line[0]
    version_starts= len(line) - line[::-1].find(' ')
    line= line[version_starts:len(line)-1]
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_ref_genome(report_list):
    tag= 'Bowtie was run against the bisulfite genome of '
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    line= line[0]
    ref_genome= re.sub(tag, '', line)
    line= re.sub(' with the specified options: .*', '', ref_genome)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_bowtie_options(report_list):
    tag= ' with the specified options: '
    line= [x for x in report_list if tag in x]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    line= line[0]
    line= line[line.find(tag) + len(tag): ]
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_is_directional(report_list):
    tag= "Option '--directional' specified: "
    line= [x for x in report_list if tag in x]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 1:
        return(re.sub('^get_', '',  inspect.stack()[0][3]), 'Yes')
    else:
        return(re.sub('^get_', '',  inspect.stack()[0][3]), 'No')

def get_seqs_tot(report_list):
    tag= 'Sequence pairs analysed in total:\t'
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        tag= 'Sequences analysed in total:\t'
        line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_seqs_alned_uniq(report_list):
    tag= 'Number of paired-end alignments with a unique best hit:\t'
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        tag= 'Number of alignments with a unique best hit from the different alignments:\t'
        line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_mapping_efficiency(report_list):
    tag= 'Mapping efficiency:\t'
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line).rstrip('%')
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_seqs_not_alned(report_list):
    tag= 'Sequence pairs with no alignments under any condition:\t'
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        tag= 'Sequences with no alignments under any condition:\t'
        line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_seqs_not_uniq(report_list):
    tag= 'Sequence pairs did not map uniquely:\t'
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        tag= 'Sequences did not map uniquely:\t'
        line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_discarded_seqs(report_list):
    tag= 'Sequence pairs which were discarded because genomic sequence could not be extracted:\t'
    line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        tag= 'Sequences which were discarded because genomic sequence could not be extracted:\t'
        line= [x for x in report_list if x.startswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_OT(report_list):
    tag= '\t((converted) top strand)'
    line= [x for x in report_list if x.endswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    line= line.split('\t')[1]
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_OB(report_list):
    tag= '\t((converted) bottom strand)'
    line= [x for x in report_list if x.endswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    line= line.split('\t')[1]
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_CTOT(report_list):
    tag= '\t(complementary to (converted) top strand)'
    line= [x for x in report_list if x.endswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    line= line.split('\t')[1]
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_CTOB(report_list):
    tag= '\t(complementary to (converted) bottom strand)'
    line= [x for x in report_list if x.endswith(tag)]
    if line == []:
        sys.exit('%s: Not matching line found!' %(inspect.stack()[0][3], ))
    if len(line) > 1:
        sys.exit('%s: More than one line found!' %(inspect.stack()[0][3], ))
    line= line[0]
    line= re.sub(tag, '', line)
    line= line.split('\t')[1]
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
        values= [get_bismark_report_for(report_list), ## Each item is tuple (header_name, value)
                get_full_path(report_list),
                get_bismark_version(report_list),
                get_ref_genome(report_list),
                get_bowtie_options(report_list),
                get_is_directional(report_list),
                get_seqs_tot(report_list),
                get_seqs_alned_uniq(report_list),
                get_mapping_efficiency(report_list),
                get_seqs_not_alned(report_list),
                get_seqs_not_uniq(report_list),
                get_discarded_seqs(report_list),
                get_OT(report_list),
                get_OB(report_list),
                get_CTOT(report_list),
                get_CTOB(report_list),
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

