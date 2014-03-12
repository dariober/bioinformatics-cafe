#!/usr/bin/env python

import subprocess
import sys
import os
import inspect
import genome_graphs
import re

docstring= """
DESCRIPTION

    Script to test and showcase usage of genomeGraphs

USAGE:
    genomeGraphsDemo.py <output-dir>

    output-dir:
        Dir for output files. Created if not exists, default to cuurent dir

genomeGraphs is expected to be on the PATH, as per installation
"""

def pdfOpen(pdf):
    try:
        p= subprocess.Popen(['xdg-open', pdf])
        p.wait()
        p.returncode
    except:
        try:
            p= subprocess.Popen(['open', pdf], stderr= subprocess.PIPE, stdout= subprocess.PIPE)
            p.wait()
            rt= p.returncode
        except:
            print("Can't open %s" %(pdf))

if len(sys.argv) != 2 or sys.argv[1] in ['-h', '--help']:
    print(docstring)
    sys.exit()

outdir= sys.argv[1]
if not os.path.exists(outdir):
    os.makedirs(outdir)
    
basedir= os.path.split(inspect.getfile(genome_graphs))[0]
datadir= os.path.join(basedir, 'demo')

# -----------------------------------------------------------------------------
# CHUNK 1
region= ['Q36', '0', '120']

ibam= ' '.join([os.path.join('$dd', 'grm022.gtf'),
                os.path.join('$dd', 'grm022.bam'),
                os.path.join('$dd', 'grm022.bed'),
                os.path.join('$dd', 'grm022.bedGraph')])

cmd="""dd="%(datadir)s"
echo '%(region)s' | genomeGraphs -b - \
    -i %(ibam)s \
    --ylab NA 'Read count' NA 'Read profile' NA \
    --mar 5 \
    --rcode NA \"rect(xleft= c(13, 72), xright= c(33, 83), ybottom= -1E9, ytop= 1E9, lwd= 0.5, border= 'darkred', col= '#7FFFD440')\" \
    -d %(outdir)s
""" %{'datadir': datadir, 'region': ' '.join(region), 'ibam': ibam, 'outdir': outdir}

outPdf= os.path.join(outdir, '_'.join(region) + '.pdf')
# cmdPrint= re.sub(datadir + '/', '', cmd)
cmdPrint= re.sub('    ', ' \ \n', cmd)

print('\n\n' + '-' * 80 + '\n')
print(cmdPrint)
print('Input files are in "%s"' %(datadir))
print('Output file "%s"' %(outPdf))
print('-' * 80)

p= subprocess.Popen(cmd, shell= True)
p.wait()
pdfOpen(outPdf)

# -----------------------------------------------------------------------------
# CHUNK 2
region= ['Q36', '60', '120']
ibam= ' '.join([os.path.join(datadir, 'grm022.gtf'),
                os.path.join(datadir, 'grm022.bam'),
                os.path.join(datadir, 'grm022.bed'),
                os.path.join(datadir, 'grm022.bedGraph')])

cmd="""echo '%(region)s' | genomeGraphs -b - \
    -i %(ibam)s \
    --ylab NA 'Read count' NA 'Read profile' NA \
    --mar 5 \
    -d %(outdir)s
""" %{'region': ' '.join(region), 'ibam': ibam, 'outdir': outdir}

outPdf= os.path.join(outdir, '_'.join(region) + '.pdf')
cmdPrint= re.sub(datadir + '/', '', cmd)
cmdPrint= re.sub('    ', ' \ \n', cmdPrint)

print('\n\n' + '-' * 80 + '\n')
print(cmdPrint)
print('Input files are in "%s"' %(datadir))
print('Output file "%s"' %(outPdf))
print('-' * 80)

p= subprocess.Popen(cmd, shell= True)
p.wait()
pdfOpen(outPdf)

# -----------------------------------------------------------------------------
# CHUNK 3
region= ['Q36', '60', '120']
ibam= ' '.join([os.path.join('$dd', 'grm022.gtf')]
    + [os.path.join('$dd', 'grm022.bedGraph')] * 50
    + [os.path.join('$dd', 'grm022.gtf')])

cmd="""dd="%(datadir)s"

cols=`printf ' grey50%%.0s' {1..50}`
vh=`printf ' 1%%.0s' {1..50}`

echo '%(region)s' | genomeGraphs -b - \
    -i %(ibam)s \
    --title 'Stacking several tracks' \
    --bg white \
    --col_track pink $cols pink \
    --col_line darkorange \
    --lwd 2 \
    -vh 5 $vh 5 \
    -d %(outdir)s
""" %{'datadir': datadir, 'region': ' '.join(region), 'ibam': ibam, 'outdir': outdir}

outPdf= os.path.join(outdir, '_'.join(region) + '.pdf')
cmdPrint= re.sub('    ', ' \\\n', cmd)

print('\n\n' + '-' * 80 + '\n')
print(cmdPrint)
print('Output file "%s"' %(outPdf))
print('-' * 80)

p= subprocess.Popen(cmd, shell= True)
p.wait()
pdfOpen(outPdf)

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# CHUNK 3
region= ['Q36', '60', '120', 'vh']
ibam= ' '.join([os.path.join(datadir, 'grm022.bedGraph'),
                os.path.join(datadir, 'grm022.bedGraph'),
                os.path.join(datadir, 'grm022.bedGraph'),
                ])

cmd="""echo '%(region)s' | genomeGraphs -b - \
    -i %(ibam)s \
    --title 'Adjusting panel size with -vh' \
    --col_track grey60 \
    --lwd 2 \
    --col_line chocolate \
    -vh 3 2 1 \
    -d %(outdir)s
""" %{'region': ' '.join(region), 'ibam': ibam, 'outdir': outdir}

outPdf= os.path.join(outdir, '_'.join(region) + '.pdf')
cmdPrint= re.sub(datadir + '/', '', cmd)
cmdPrint= re.sub('    ', ' \ \n', cmdPrint)

print('\n\n' + '-' * 80 + '\n')
print(cmdPrint)
print('Input files are in "%s"' %(datadir))
print('Output file "%s"' %(outPdf))
print('-' * 80)

p= subprocess.Popen(cmd, shell= True)
p.wait()
pdfOpen(outPdf)

# -----------------------------------------------------------------------------

sys.exit()

"""
dd="/Users/berald01/svn_checkout/bioinformatics-misc/branches/genomeGraphs/genome_graphs/demo"

bdg=`printf " $dd/grm022.bedGraph%.0s" {1..100}`
cols=`printf ' red blue%.0s' {1..50}`
vh=`printf ' 1%.0s' {1..100}`

echo 'Q36 0 1000' | genomeGraphs -b - \
 -i $dd/grm022.gtf $bdg $dd/grm022.gtf \
 --title 'Stacking several tracks' \
 --bg white \
 --col_track pink $cols pink \
 --col_line NA \
 --col_mark NA \
 --names '' \
 --col_grid NA \
 --fbg white \
 --lwd 2 \
 -vh 4 $vh 4 \
 -mh 3 10 \
 -W 16 \
 -H 24 \
 --col_yaxis NA \
 --cex_axis 1.5 \
 -d tmp
"""