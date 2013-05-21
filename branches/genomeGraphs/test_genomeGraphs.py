#!/usr/bin/env py.test

"""Run me:  
py.test ~/svn_checkout/bioinformatics-misc/branches/genomeGraphs/test_genomeGraphs.py
"""
import shutil
import os
import sys
import subprocess as sp

genomeGraphs='~/svn_checkout/bioinformatics-misc/branches/genomeGraphs/genomeGraphs.py'
example_dir='/Users/berald01/Documents/genomeGraphs/example'
os.chdir(example_dir)
tmpdir= 'wdir'
outdir= 'pdf'

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
    for x in ['depth', 'A', 'a', 'C', 'c', 'G', 'g', 'T', 't', 'N', 'n', 'Z', 'z']:
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
    assert line[n*1] == '922'      ## Depth
    assert line[n*2] == '0'       ## A
    assert line[n*3] == '0'       ## a
    assert line[n*4] == '510'       ## C
    assert line[n*5] == '411'       ## c
    assert line[n*6] == '0'       ## G
    assert line[n*7] == '0'       ## g
    assert line[n*8] == '0'      ## T
    assert line[n*9] == '1'       ## t
    assert line[n*10] == '0'       ## N
    assert line[n*11] == '0'       ## n
    assert line[n*12] == '510'      ## Z: Sum ACGTN.
    assert line[n*13] == '412'      ## z: Sum actgn.

    line= pileup[556]
    print(line)
    assert line[2] == '5567333' ## Make sure you are at this position

    ## Second bam
    assert line[(n*1)+1] == '462'      ## Depth
    assert line[(n*2)+1] == '0'       ## A
    assert line[(n*3)+1] == '0'       ## a
    assert line[(n*4)+1] == '0'       ## C
    assert line[(n*5)+1] == '1'       ## c
    assert line[(n*6)+1] == '0'       ## G
    assert line[(n*7)+1] == '0'       ## g
    assert line[(n*8)+1] == '227'      ## T
    assert line[(n*9)+1] == '234'       ## t
    assert line[(n*10)+1] == '0'       ## N
    assert line[(n*11)+1] == '0'       ## n
    assert line[(n*12)+1] == '227'      ## Z: Sum ACGTN.
    assert line[(n*13)+1] == '235'      ## z: Sum actgn.
#    ## Third bam
#    assert line[(n*1) + 2] == '12'
#    assert line[(n*2) + 2] == '0'
#    assert line[(n*3) + 2] == '2'
#    assert line[(n*4) + 2] == '0'
#    assert line[(n*5) + 2] == '10'
#    assert line[(n*6) + 2] == '0'
#    assert line[(n*7) + 2] == '12'
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
    ## _MEMO_: On shell remember to use `echo -e` to correctly read \t!!
    cmd= """echo 'chr7\t5567130\t5567189' | %(genomeGraphs)s -i bam/ds051.actb.bam -b - -f annotation/chr7.fa --tmpdir %(tmpdir)s -d %(outdir)s"""  %{'genomeGraphs': genomeGraphs, 'tmpdir': tmpdir, 'outdir':outdir}
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    seqbed= open(os.path.join(tmpdir, 'chr7_5567130_5567189.seq.txt')).readlines()
    seqbed= [x.strip().split('\t') for x in seqbed]
    assert seqbed[0] == ['chrom', 'start', 'end', 'base'] ## Check first line is this header
    nucseq= ''.join([x[3] for x in seqbed[1:]]) ## Get column 'base', skipping header.
    assert nucseq == 'TGACTATTAAAAAAACAACAATGTGCAATCAAAGTCCTCGGCCACATTGTGAACTTTGGG' ## This sequence manually checked.
