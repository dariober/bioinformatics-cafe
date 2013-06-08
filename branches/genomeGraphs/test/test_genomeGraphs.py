#!/usr/bin/env py.test

"""Run me:  
py.test ~/svn_checkout/bioinformatics-misc/branches/genomeGraphs/test_genomeGraphs.py
"""
import shutil
import os
import sys
import subprocess as sp
import pybedtools
import pycoverage
import inspect

genomeGraphs='~/svn_checkout/bioinformatics-misc/branches/genomeGraphs/genomeGraphs.py'
example_dir='/Users/berald01/Documents/genomeGraphs/example'
os.chdir(example_dir)
tmpdir= 'pytest_wdir'
outdir= 'pytest_pdf'

def test_only_one_gtf():
    cmd= '%(genomeGraphs)s -i annotation/genes.gtf.gz -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    print(cmd)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_only_one_bedgraph():
    cmd= '%(genomeGraphs)s -i bedgraph/profile.bedGraph.gz -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    print(cmd)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_only_one_bam():
    cmd= '%(genomeGraphs)s -i bam/ds051.actb.bam -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    print(cmd)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_all_bams():
    cmd= '%(genomeGraphs)s -i bam/ds*.bam -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    print(cmd)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_bam_gtf_bedgraph():
    cmd= '%(genomeGraphs)s -i annotation/genes.gtf.gz bam/ds*.bam  annotation/genes.gtf.gz bedgraph/profile.bedGraph.gz annotation/genes.gtf.gz -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_file_headers():
    """Make sure headers are what you think they are
    """
    cmd= 'echo "chr7\t5566757\t5566829" | %(genomeGraphs)s -i annotation/genes.gtf.gz bam/ds*.bam  annotation/genes.gtf.gz bedgraph/profile.bedGraph.gz annotation/genes.gtf.gz -b - -f annotation/chr7.fa --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    print(cmd)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == '' 
    header= open(os.path.join(tmpdir, 'chr7_5566757_5566829.seq.txt')).readline().strip().split('\t')
    assert header == ['chrom', 'start', 'end', 'base']

    header= open(os.path.join(tmpdir, 'chr7_5566757_5566829.grp.bed.txt')).readline().strip().split('\t')

    expected_header= ['chrom', 'start', 'end']
    for x in ['A', 'a', 'C', 'c', 'G', 'g', 'T', 't', 'N', 'n', 'Z', 'z']:
        for bam in ['bam/ds051.actb.bam', 'bam/ds052.actb.bam', 'bam/ds053.actb.bam']:
            expected_header.append(bam + '.' + x)
    assert header == expected_header

def test_bam_coverage_eq_igv():
    """Confirm that the counts obtained from the BAM files is consistent with
    IGV.
    This test is very important: It means the counts in *.grp.bed.txt are correct.
    """
    n= 3 ## No. BAMs you have in input
    cmd= '%(genomeGraphs)s -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam -b actb.bed --nwinds 5000 --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    print(cmd)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    pileup= open(os.path.join(tmpdir, 'chr7_5566755_5567571_ACTB.grp.bed.txt')).readlines()
    pileup= [x.strip().split('\t') for x in pileup]
    line= pileup[599]
    print(line)
    assert line[2] == '5567376' ## Make sure you are at this position
    """These counts cross-checked in IGV.
    Multiply by the number of bam files you have in `-i` to get the counts
    """
#    assert line[n*1] == '922'      ## Depth
    assert line[n*1] == '0'       ## A
    assert line[n*2] == '0'       ## a
    assert line[n*3] == '510'       ## C
    assert line[n*4] == '411'       ## c
    assert line[n*5] == '0'       ## G
    assert line[n*6] == '0'       ## g
    assert line[n*7] == '0'      ## T
    assert line[n*8] == '1'       ## t
    assert line[n*9] == '0'       ## N
    assert line[n*10] == '0'       ## n
    assert line[n*11] == '510'      ## Z: Sum ACGTN.
    assert line[n*12] == '412'      ## z: Sum actgn.

    line= pileup[556]
    print(line)
    assert line[2] == '5567333' ## Make sure you are at this position

    ## Second bam
#    assert line[(n*1)+1] == '462'      ## Depth
    assert line[(n*1)+1] == '0'       ## A
    assert line[(n*2)+1] == '0'       ## a
    assert line[(n*3)+1] == '0'       ## C
    assert line[(n*4)+1] == '1'       ## c
    assert line[(n*5)+1] == '0'       ## G
    assert line[(n*6)+1] == '0'       ## g
    assert line[(n*7)+1] == '227'      ## T
    assert line[(n*8)+1] == '234'       ## t
    assert line[(n*9)+1] == '0'       ## N
    assert line[(n*10)+1] == '0'       ## n
    assert line[(n*11)+1] == '227'      ## Z: Sum ACGTN.
    assert line[(n*12)+1] == '235'      ## z: Sum actgn.

    if stderr == '':
        shutil.rmtree(tmpdir)

def test_annotation_bed():
    """Confirms a bed file intersect at the right positions
    """
    cmd= '%(genomeGraphs)s -i annotation/actb_ann.bed -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    print(cmd)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    annbed= open(os.path.join(tmpdir, 'chr7_5566755_5567571_ACTB.nonbam.bed.txt')).readlines()
    annbed= [x.strip().split('\t') for x in annbed]
    assert len(annbed) == 4
    
    igv_starts= ['5566844', '5567044', '5567244', '5567444'] ## This coords extracted from IGV 
    igv_end= ['5566944', '5567144', '5567344', '5567544']    ## after loading actb_ann.bed.
    for i in range(0, len(annbed)):
        line= annbed[i]
        start= line[1]
        end= line[2]
        assert start ==  igv_starts[i]
        assert end ==  igv_end[i]
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_annotation_gtf():
    """Confirms a bed file intersect at the right positions
    """
    cmd= '%(genomeGraphs)s -i annotation/genes.gtf.gz -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    annbed= open(os.path.join(tmpdir, 'chr7_5566755_5567571_ACTB.nonbam.bed.txt')).readlines()
    annbed= [x.strip().split('\t') for x in annbed]
    assert len(annbed) == 3

    gtf_starts= ['5566777', '5567377', '5567380']
    gtf_end= ['5567522', '5567381', '5567522'] 
    for i in range(0, len(annbed)):
        line= annbed[i]
        start= line[1]
        end= line[2]
        assert start ==  gtf_starts[i]
        assert end ==  gtf_end[i]
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_plot_params():
    """Test a number of plotting parameters.
    """
    cmd= '%(genomeGraphs)s -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

    """Use --replot: change size
    """

    cmd= '%(genomeGraphs)s --replot -W 15 -H 18 -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

    """Use --replot: change proportions
    """
    cmd= '%(genomeGraphs)s --replot --vheights 1 2 2 2 3 1 -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

    """Use --replot: change labels sizes
    """
    cmd= '%(genomeGraphs)s --replot --cex 2 -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

    """Use --replot:
    - Assign plot names other than default file names. Put one name less: It should be recycled.
    - Assign colour to names (recycled)
    - Change size 
    """
    cmd= '%(genomeGraphs)s --replot --names gtf ds051 ds052 ds053 bedgraph --col_names red blue --cex_names 1.75 -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

    """Change track colour
    """
    cmd= '%(genomeGraphs)s --replot --col_track red blue blue green -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %(tmpdir)s -d %(outdir)s' %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

def test_sequence_and_pipe():
    cmd= """echo 'chr7\t5567130\t5567189' | %(genomeGraphs)s -i bam/ds051.actb.bam -b - --slop 0 0 -f annotation/chr7.fa --tmpdir %(tmpdir)s -d %(outdir)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    seqbed= open(os.path.join(tmpdir, 'chr7_5567130_5567189.seq.txt')).readlines()
    seqbed= [x.strip().split('\t') for x in seqbed]
    assert seqbed[0] == ['chrom', 'start', 'end', 'base'] ## Check first line is this header
    nucseq= ''.join([x[3] for x in seqbed[1:]]) ## Get column 'base', skipping header.
    print(nucseq)
    assert nucseq == 'GACTATTAAAAAAACAACAATGTGCAATCAAAGTCCTCGGCCACATTGTGAACTTTGGG' ## This sequence manually checked.

def test_sequence_and_pipe2():
    cmd= """echo 'chr7 5567130 5567189' | %(genomeGraphs)s -i bam/ds051.actb.bam -b - --slop 1 1 -f annotation/chr7.fa --tmpdir %(tmpdir)s -d %(outdir)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    seqbed= open(os.path.join(tmpdir, 'chr7_5567130_5567189.seq.txt')).readlines()
    seqbed= [x.strip().split('\t') for x in seqbed]
    assert seqbed[0] == ['chrom', 'start', 'end', 'base'] ## Check first line is this header
    nucseq= ''.join([x[3] for x in seqbed[1:]]) ## Get column 'base', skipping header.
    assert nucseq == 'TGACTATTAAAAAAACAACAATGTGCAATCAAAGTCCTCGGCCACATTGTGAACTTTGGGG' ## This sequence manually checked.


def test_slopbed_pybedtool():
    """Test function slopbed with pybedtool interval feature 
    """
    f = iter(pybedtools.BedTool('chr1 1 100 asdf 0 + a b c d', from_string=True)).next()
    xf= pycoverage.slopbed(f, [0, 0])
    assert f == xf
    
    xf= pycoverage.slopbed(f, [0, 10])
    assert xf == iter(pybedtools.BedTool('chr1 1 110 asdf 0 + a b c d', from_string=True)).next()
    
    xf= pycoverage.slopbed(f, [1, 0])
    assert xf == iter(pybedtools.BedTool('chr1 0 100 asdf 0 + a b c d', from_string=True)).next()

    xf= pycoverage.slopbed(f, [10, 0]) ## Left coord must remain 0
    assert xf == iter(pybedtools.BedTool('chr1 0 100 asdf 0 + a b c d', from_string=True)).next()

    xf= pycoverage.slopbed(f, [1, 1]) ## Left coord must remain 0
    assert xf == iter(pybedtools.BedTool('chr1 0 101 asdf 0 + a b c d', from_string=True)).next()

    xf= pycoverage.slopbed(f, [0.1, 0.13]) ## Left coord must remain 0
    assert xf == iter(pybedtools.BedTool('chr1 0 113 asdf 0 + a b c d', from_string=True)).next()

    xf= pycoverage.slopbed(f, [0.1, 1.13]) ## Left coord must remain 0
    assert xf == iter(pybedtools.BedTool('chr1 0 212 asdf 0 + a b c d', from_string=True)).next()

    try:
        xf= pycoverage.slopbed(f, [-1, 0])
    except pycoverage.SlopError:
        assert True
    try:
        xf= pycoverage.slopbed(f, [0, '0'])
    except pycoverage.SlopError:
        assert True
        
def test_slopbed_list():
    """Test function slopbed with list interval feature 
    """
    f= ['chr1', 1, 100, 'actb']
    xf= pycoverage.slopbed(f, [0, 0])
    assert f == xf

    xf= pycoverage.slopbed(f, [0.1, 1.13]) ## Left coord must remain 0
    assert xf == ['chr1', 0, 212, 'actb']

    try:  
        xf= pycoverage.slopbed(['chr1', '0', 212, 'actb'], [0.1, 1.13]) ## Left coord must remain 0
    except TypeError:
        assert True

def test_stdin_inbed_to_fh():
    inbed= open('actb.bed')
    bedfh= pycoverage.stdin_inbed_to_fh(inbed)
    inlist= open(bedfh.name).readlines()
    os.remove(bedfh.name)
    assert len(inlist) == 5
    assert inlist[1] == 'chr7\t5500000\t5566830\tACTB\n'

def test_overplotting_RGB_in_one_graph():
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        -op 1 1 1 \
        --col_line red green blue \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0

def test_overplotting_RB_in_one_graph_G_alone():
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        --title 'Red and blue in top plot. Green alone. R at ~1.5; B at ~3.5; G at ~2.5' \
          -op        1   2     1 \
        --col_line red green blue \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    print(stderr)
    assert p.returncode == 0

def test_overplotting_G_alone_RB_in_one_graph():
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        -op        2   1     2 \
        --col_line red green blue \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0

def test_overplotting_G_alone_RB_in_one_graph_ymax():
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        --ymax max \
        -op        2   1     2 \
        --col_line red green blue \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0

def test_overplotting_G_alone_RB_in_one_graph_ylim_to_5_and_4():
    """Memo: the plotting order is given by -op not by the input list.
    """
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        --ymax         4     5 \
        -op        2   1     2 \
        --col_line red green blue \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0

def test_overplotting_RGB_in_one_graph_annot():
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        -op 1 1 1 \
        --col_line red green blue \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
           actb.bed \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0

def test_vheight_all_equal():
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        -vh 1 1 1 1 \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
           actb.bed \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0

def test_vheight_4x4x2x1x():
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        -vh 4 4 2 1 \
        --title 'First and second plot 2x the third and 4x the fourth' \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
           actb.bed \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0

def test_vheight_4x1x_recycled():
    """1st 4x the 2nd and 3rd 4x the 4th.
    """
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        -vh 4 1 \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
           actb.bed \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0
    
def test_mar_height_1top_1bottom():
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        -mh 1 1 \
        --title ':region: <- region coords here\nTop and bottom margin the same height.\nQuite large.' \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
           actb.bed \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0

def test_mar_height_01top_03bottom_cex_resized():
    outfile= inspect.stack()[0][3] + '.pdf'
    cmd= """%(genomeGraphs)s \
        -mh 0.2 0.6 \
        --title 'Bottom mar 3x top.\nTop cex resized because on two lines' \
        -i bedgraph/profile.odd.bedGraph \
           bedgraph/profile.odd_plus1.bedGraph \
           bedgraph/profile.odd_plus2.bedGraph \
           actb.bed \
        -b actb.bed --tmpdir %(tmpdir)s -o %(outfile)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outfile': os.path.join(outdir, outfile)}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert p.returncode == 0