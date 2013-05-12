#!/usr/bin/env test.py

"""Run me:
py.test ~/svn_checkout/bioinformatics-misc/trunk/genomeGraphs/test_genomeGraphs.py
"""
import shutil
import os
import sys
import subprocess as sp

example_dir='/Users/berald01/Documents/genomeGraphs/example'
os.chdir(example_dir)
tmpdir= 'wdir'
outdir= 'pdf'

def test_only_one_gtf():
    cmd= 'genomeGraphs.py -i annotation/genes.gtf.gz -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_only_one_bedgraph():
    cmd= 'genomeGraphs.py -i bedgraph/profile.bedGraph.gz -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_only_one_bam():
    cmd= 'genomeGraphs.py -i bam/ds051.actb.bam -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_all_bams():
    cmd= 'genomeGraphs.py -i bam/*.bam  -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_bam_gtf_bedgraph():
    cmd= 'genomeGraphs.py -i annotation/genes.gtf.gz bam/*.bam  annotation/genes.gtf.gz bedgraph/profile.bedGraph.gz annotation/genes.gtf.gz -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_bam_coverage_eq_igv():
    """Confirm that the counts obtained from the BAM files is consistent with
    IGV.
    """
    n= 3 ## No. BAMs you have in input
    cmd= 'genomeGraphs.py -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''
    pileup= open(os.path.join(tmpdir, 'chr7_5566757_5566829_ACTB.grp.bed.txt')).readlines()
    pileup= [x.strip().split('\t') for x in pileup]
    line= pileup[16]
    print(line)
    assert line[2] == '5566793' ## Make sure you are at this position
    """These counts cross-checked in IGV.
    Multiply by the number of bam files you have in `-i` to get the counts
    """
    assert line[n*1] == '22'      ## Depth
    assert line[n*2] == '0'       ## A
    assert line[n*3] == '1'       ## C
    assert line[n*4] == '1'       ## G
    assert line[n*5] == '20'      ## T
    assert line[n*6] == '0'       ## N
    assert line[n*7] == '22'      ## Sum ACGTN.
    ## Second bam
    assert line[(n*1) + 1] == '12'
    assert line[(n*2) + 1] == '0'
    assert line[(n*3) + 1] == '6'
    assert line[(n*4) + 1] == '0'
    assert line[(n*5) + 1] == '6'
    assert line[(n*6) + 1] == '0'
    assert line[(n*7) + 1] == '12'
    ## Third bam
    assert line[(n*1) + 2] == '12'
    assert line[(n*2) + 2] == '0'
    assert line[(n*3) + 2] == '2'
    assert line[(n*4) + 2] == '0'
    assert line[(n*5) + 2] == '10'
    assert line[(n*6) + 2] == '0'
    assert line[(n*7) + 2] == '12'
    if stderr == '':
        shutil.rmtree(tmpdir)

def test_annotation_bed():
    """Confirms a bed file intersect at the right positions
    """
    cmd= 'genomeGraphs.py -i annotation/actb_ann.bed -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
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
    cmd= 'genomeGraphs.py -i annotation/genes.gtf.gz -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
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
    cmd= 'genomeGraphs.py -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

    """ Use --replot: change size
    """
    cmd= 'genomeGraphs.py --replot -W 15 -H 18 -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

    """Use --replot: change proportions
    """
    cmd= 'genomeGraphs.py --replot --vheights 1 2 2 2 3 1 -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

    """ Use --replot: change labels sizes
    """
    cmd= 'genomeGraphs.py --replot --cex 2 -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

    """ Use --replot:
    - Assign plot names other than default file names.
    Put one name less: It should be recycled.
    - Assign colour to names (recycled)
    - Change size 
    """
    cmd= 'genomeGraphs.py --replot --names gtf ds051 ds052 ds053 bedgraph --col_names red blue --cex_names 1.75 -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

    """Change track colour
    """
    cmd= 'genomeGraphs.py --replot --col_track red blue blue green -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    assert stderr == ''

#def test_plot_params():
#    """Test a number of plotting parameters.
#    """
#    cmd= 'genomeGraphs.py -i annotation/genes.gtf.gz bam/ds051.actb.bam bam/ds052.actb.bam bam/ds053.actb.bam bedgraph/profile.bedGraph.gz annotation/actb_ann.bed -b actb.bed --tmpdir %s -d %s' %(tmpdir, outdir)
#    p= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
#    stdout, stderr= p.communicate()
#    assert stderr == ''
