#!/usr/bin/python

import re
import sys

psq_re_f= re.compile('G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}')
psq_re_r= re.compile('C{3,}.{1,7}C{3,}.{1,7}C{3,}.{1,7}C{3,}')

ref_seq_fh= open(sys.argv[1])

ref_seq=[]
line= (ref_seq_fh.readline()).strip()
chr= re.sub('^>', '', line)
line= (ref_seq_fh.readline()).strip()
while True:
    while line.startswith('>') is False:
        ref_seq.append(line)
        line= (ref_seq_fh.readline()).strip()
        if line == '':
            break
    ref_seq= ''.join(ref_seq)
    for m in re.finditer(psq_re_f, ref_seq):
        quad_id= str(chr) + '_' + str(m.start()) + '_' + str(m.end()) + '_for'
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s' %(chr, m.start(), m.end(), quad_id, len(m.group(0)), '+', m.group(0)))
    for m in re.finditer(psq_re_r, ref_seq):
        quad_id= str(chr) + '_' + str(m.start()) + '_' + str(m.end()) + '_rev'
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s' %(chr, m.start(), m.end(), quad_id, len(m.group(0)), '+', m.group(0)))
    chr= re.sub('^>', '', line)
    ref_seq= []
    line= (ref_seq_fh.readline()).strip()
    if line == '':
        break
sys.exit()
