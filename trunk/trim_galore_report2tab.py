#!/home/berald01/.local/bin/python

import argparse
import sys
import os
import re
import inspect
import traceback
import subprocess

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Parse one or more reports from trim_galore to tabular format.
    Formatted output can be uplaoded to sblab. Only rows not already found in trim_galore_report
    are inserted.
EXAMPLE
    ## Parse all reports in pwd:
    ls *_trimming_report.txt | trim_galore_report2tab.py --first_header -
    ## Upload to sblab
    trim_galore_report2tab.py sh023_hela_i_trimming_report.txt --sblab_upload
TODO:
    - For uplaoding to sblab: Name of data file and sql script file are hardcoded.
      This means that if multiple instances of trim_galore_report2tab are executed
      simulataneously  the data and ssql files might be overwritten.
      
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
parser.add_argument('--sblab_upload',
                   action= 'store_true',
                   help='''Uplaod the report(s) to sblab database.
With this option no output to stdin other than the upload logs. It is not necessary to
to be on the machine hosting sblab, but the local machine must be abple to ssh/scp to
sblab machine
                   ''')

args = parser.parse_args()

# -----------------------------------------------------------------------------
def list2pgarray(tuplelist):
    """Convert a list of tuple to a string formatted like an 2D-array suiatble for postgres
    list2pgarray([('3', '565296'), ('4', '113096'), (5, 30513.0)]) => "{{'3', '565296'}, {'4', '113096'}, {5, 30513.0}}"
    The tuples are expected to be 
    """
    strarray= str(tuplelist).replace(')', '}').replace('(', '{').replace('[', '{').replace(']', '}')
    return(strarray)    
# -----------------------------------------------------------------------------

def get_input_filename(report_list):
    tag= 'Input filename: '
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_phred_cutoff(report_list):
    tag= 'Quality Phred score cutoff: '
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_phred_cutoff(report_list):
    tag= 'Quality encoding type selected: '
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= re.sub(tag, '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_adapter_seq(report_list):
    tag= "Adapter sequence: "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= re.sub(tag, '', line).strip("'")
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_min_adapt_overlap(report_list):
    tag= "Minimum required adapter overlap (stringency): "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '').rstrip(" bp")
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_is_paired(report_list):
    tag= 'Minimum required sequence length for both reads before a sequence pair gets removed: '
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) == 1:
        return(re.sub('^get_', '',  inspect.stack()[0][3]), True)
    tag= 'Minimum required sequence length before a sequence gets removed: '
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) == 1:
        return(re.sub('^get_', '',  inspect.stack()[0][3]), False)
    sys.exit('%s: Cannot determine whether file is paired' %(inspect.stack()[0][3], ))

def get_cutadapt_version(report_list):
    tag= "cutadapt version "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '')
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_cutadapt_parameters(report_list):
    tag= "Command line parameters: "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '')
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_max_error_perc(report_list):
    tag= "Maximum error rate: "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '').strip('%')
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_proc_reads(report_list):
    tag= "Processed reads: "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '').strip('%')
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_trimmed(report_list):
    tag= "Trimmed reads: "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '')
    line= re.sub(' .*', '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_trimmed_perc(report_list):
    tag= "Trimmed reads: "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '').rstrip('%)')
    line= re.sub('\d+ \( ', '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_too_short(report_list):
    tag= "Too short reads: "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '')
    line= re.sub(' .*', '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_too_short_perc(report_list):
    tag= "Too short reads: "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '').rstrip('% of processed reads)')
    line= re.sub('\d+ \( ', '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_time_sec(report_list):
    tag= "Total time: "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '').rstrip(' s')
    line= re.sub(' *', '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_time_read_ms(report_list):
    tag= "Time per read: "
    line= [x for x in report_list if x.startswith(tag)]
    if len(line) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    if len(line) == 0:
        sys.exit('%s: No line found for "%s"' %(inspect.stack()[0][3], tag))
    line= line[0]
    line= line.replace(tag, '').rstrip(' ms')
    line= re.sub(' *', '', line)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), line)

def get_adapt_hist(report_list):
    tag= 'Histogram of adapter lengths'
    hist_start= report_list.index(tag) + 2
    hist_end= [x for x in report_list if x.startswith('RUN STATISTICS FOR INPUT FILE: ')]
    if len(hist_end) > 1:
        sys.exit('%s: More than one line found' %(inspect.stack()[0][3], ))
    hist_end= report_list.index(hist_end[0])
    hist= report_list[hist_start:hist_end]
    hist2= [tuple(x.split('\t')) for x in hist]
    hist= []
    for t in hist2:
        t= (int(t[0]), int(t[1]))
        hist.append(t)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), hist)

def get_adapt_hist_freq(report_list):
    hist= get_adapt_hist(report_list)[1]
    counts= [int(x[1]) for x in hist]
    lengths= [x[0] for x in hist]
    nreads= float(get_proc_reads(report_list)[1])
    hist_freq= [round((x/nreads)*100, 2) for x in counts]
    hist_freq2= zip(lengths, hist_freq)
    hist_freq= []
    for t in hist_freq2:
        t= (int(t[0]), float(t[1]))
        hist_freq.append(t)
    return(re.sub('^get_', '',  inspect.stack()[0][3]), hist_freq)
# -----------------------------------------------------------------------------

if args.infile == ['-']:
    args.infile= sys.stdin.readlines()
    args.infile= [x.strip() for x in args.infile]

data= [] ## Store each report as string here as list
for trim_report in args.infile:
    report_list= open(trim_report).readlines()
    report_list= [x.strip() for x in report_list if not x.strip() == '']
    try:
        adapt_hist= get_adapt_hist(report_list)
        adapt_hist= (adapt_hist[0], list2pgarray(adapt_hist[1]))
        adapt_hist_freq= get_adapt_hist_freq(report_list)
        adapt_hist_freq= (adapt_hist_freq[0], list2pgarray(adapt_hist_freq[1]))
        values= [get_input_filename(report_list), ## Each item is tuple (header_name, value)
                get_phred_cutoff(report_list),
                get_adapter_seq(report_list),
                get_min_adapt_overlap(report_list),
                get_is_paired(report_list),
                get_cutadapt_version(report_list),
                get_cutadapt_parameters(report_list),
                get_max_error_perc(report_list),
                get_proc_reads(report_list),
                get_trimmed(report_list),
                get_trimmed_perc(report_list),
                get_too_short(report_list),
                get_too_short_perc(report_list),
                get_time_sec(report_list),
                get_time_read_ms(report_list),
                adapt_hist,
                adapt_hist_freq,
                ('report_name', trim_report)]
    except:
        stack_trace = traceback.format_exc()
        print('\nReport "%s" raised an exception:\n' %(trim_report))
        sys.exit(stack_trace)
        
    data.append(values)

if args.sblab_upload:
    reportname= 'trim_galore_report_upload.tsv'
    fout= open(reportname, 'w')
    for report in data:
        fout.write('\t'.join([str(v[1]) for v in report]) + '\n')
    fout.close()
    
    sql= """-- Script to upload trim_galore_report: 
CREATE TEMP TABLE tmp_trim_galore_report (LIKE trim_galore_report);
COPY tmp_trim_galore_report FROM '%s' WITH CSV DELIMITER E'\\t' NULL AS 'NA';
INSERT INTO trim_galore_report (SELECT * FROM tmp_trim_galore_report EXCEPT SELECT * FROM trim_galore_report);
DROP TABLE tmp_trim_galore_report; """ %("/tmp/" + reportname)
    sqlscript= reportname + '.upload.tmp.sql'
    fout= open(sqlscript, 'w')
    fout.write(sql)
    fout.close()
    ## Send data file and SQL script to sblab machine
    cmd= 'scp %s $mac_office:/tmp; scp %s $mac_office:/tmp' %(reportname, sqlscript)
    print(cmd)
    p= subprocess.Popen(cmd, shell= True)
    p.wait()
    
    ## Upload data
    cmd= '''ssh $mac_office 'source ~/.bash_profile;
    chmod 777 /tmp/%(report)s;
    psql -U dberaldi -d sblab -w < /tmp/%(sqlscript)s;
    rm /tmp/%(report)s;
    rm /tmp/%(sqlscript)s;
    exit' ''' %{'report': reportname, 'sqlscript': sqlscript}
    print(cmd)
    p= subprocess.Popen(cmd, shell= True)
    p.wait()

    os.remove(reportname)
    os.remove(sqlscript)
    sys.exit()


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