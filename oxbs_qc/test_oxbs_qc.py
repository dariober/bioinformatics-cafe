#!/usr/bin/env py.test

import oxbs_qc
import sys
import os
import subprocess
import oxbs_qc_func
import shutil
import gzip

def test_get_wdir():
    passed= False
    try:
        wdir= oxbs_qc_func.get_wdir('')
    except oxbs_qc_func.GetWdirException:
            passed= True
    assert passed
    
    passed= False
    try:
        wdir= oxbs_qc_func.get_wdir(['foo'])
    except oxbs_qc_func.GetWdirException:
            passed= True
    assert passed

    wdir= oxbs_qc_func.get_wdir('.')
    assert wdir == '.'
    wdir= oxbs_qc_func.get_wdir('test_out/test_wdir')
    assert wdir == 'test_out/test_wdir'
    assert os.path.exists(wdir) == True
    os.rmdir(wdir) ## _MEMO_: os.rmdir removes dir only if it is empty

try:
    shutil.rmtree('test_out')
except OSError:
    pass    
   
def test_task_trim_galore_SE():
    passed= False
    try:
        tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz'], outFastq= ['bla/fq1.fq.gz'], tmpdir= 'test_out', opts= '-o test_out/trim_galore', path= '')
    except oxbs_qc_func.TrimGaloreException:
        "outFastq has path attached"
        passed= True
    assert passed

def test_task_trim_galore_SE_1():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz'], outFastq= ['fq1.fq.gz'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_1', path= '')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_1/fq1.fq.gz'
    assert os.path.exists(tg['fastq'][0])
    assert os.path.exists('test_out/trim_galore_1/mjb042_oxBS_R1.fastq.gz_trimming_report.txt')

def test_task_trim_galore_SE_2():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.txt'], outFastq= ['fq1.txt'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_2', path= '')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_2/fq1.txt'
    assert os.path.exists(tg['fastq'][0])

def test_task_trim_galore_SE_3():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz'], outFastq= ['fq1.txt'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_3 --gzip', path= '')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_3/fq1.txt.gz'
    assert os.path.exists(tg['fastq'][0])

def test_task_trim_galore_SE_4():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz'], outFastq= ['fq1.txt'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_4 --dont_gzip', path= '')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_4/fq1.txt'
    assert os.path.exists(tg['fastq'][0])

def test_task_trim_galore_PE_exceptions():
    ## 2 inputs 1 output
    passed= False
    try:
        tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz', 'test_data/mjb042_oxBS_R2.fastq.gz'], outFastq= ['fq1.fq.gz'], tmpdir= 'test_out', opts= '', path= '')
    except oxbs_qc_func.TrimGaloreException:
        passed= True
    assert passed

    ## 1 inputs 2 output
    passed= False
    try:
        tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz'], outFastq= ['fq1.fq.gz', 'fq2.fq.gz'], tmpdir= 'test_out', opts= '', path= '')
    except oxbs_qc_func.TrimGaloreException:
        passed= True
    assert passed

    ## same input twice
    passed= False
    try:
        tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz', 'test_data/mjb042_oxBS_R1.fastq.gz'], outFastq= ['fq1.fq.gz', 'fq2.fq.gz'], tmpdir= 'test_out', opts= '', path= '')
    except oxbs_qc_func.TrimGaloreException:
        passed= True
    assert passed

    ## same Output twice
    passed= False
    try:
        tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz', 'test_data/mjb042_oxBS_R2.fastq.gz'], outFastq= ['fq1.fq.gz', 'fq1.fq.gz'], tmpdir= 'test_out', opts= '', path= '')
    except oxbs_qc_func.TrimGaloreException:
        passed= True
    assert passed

def test_task_trim_galore_PE_1():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz', 'test_data/mjb042_oxBS_R2.fastq.gz'], outFastq= ['fq1.fq.gz', 'fq2.fq.gz'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_pe1', path= '')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_pe1/fq1.fq.gz' ## Files are named as you expect
    assert tg['fastq'][1] == 'test_out/trim_galore_pe1/fq2.fq.gz'
    assert os.path.exists(tg['fastq'][0]) ## Files do exist
    assert os.path.exists(tg['fastq'][1])

def test_task_trim_galore_PE_2():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz', 'test_data/mjb042_oxBS_R2.txt.gz'], outFastq= ['fq1.fq.gz', 'fq2.fq.gz'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_pe2', path= '')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_pe2/fq1.fq.gz' ## Files are named as you expect
    assert tg['fastq'][1] == 'test_out/trim_galore_pe2/fq2.fq.gz'
    for x in tg['fastq']:
        assert os.path.exists(x) ## files do exist
    for x in tg['report']:
        assert os.path.exists(x) ## Reports do exist
    assert len(os.listdir('test_out/trim_galore_pe2/')) == 4

def test_task_trim_galore_PE_3():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.txt', 'test_data/mjb042_oxBS_R2.txt'], outFastq= ['fq1.txt', 'fq2.txt'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_pe3', path= '')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_pe3/fq1.txt' ## Files are named as you expect
    assert tg['fastq'][1] == 'test_out/trim_galore_pe3/fq2.txt'
    for x in tg['fastq']:
        assert os.path.exists(x) ## files do exist
    for x in tg['report']:
        assert os.path.exists(x) ## Reports do exist
    assert len(os.listdir('test_out/trim_galore_pe3/')) == 4
    
def test_trim_fastq_exception():
    passed= False
    try:
        oxbs_qc_func.trim_fastq('test_data/mjb042_oxBS_R1.txt', -1, 'mjb042_oxBS_R1.fq') ### Invalid -1
    except oxbs_qc_func.TrimFastqException:
        passed= True
    assert passed

    passed= False
    try:
        x= oxbs_qc_func.trim_fastq('test_data/mjb042_oxBS_R1.txt.gz', 1, 'test_data/mjb042_oxBS_R1.txt.gz') ## input == output
    except oxbs_qc_func.TrimFastqException:
        passed= True
    assert passed

    passed= False
    try:
        x= oxbs_qc_func.trim_fastq('test_data/mjb042_oxBS_R1.txt.gz', 1, '/nonsense/mjb042_oxBS_R1.txt.gz') ## output dir does not exists
    except oxbs_qc_func.TrimFastqException:
        passed= True
    assert passed

def test_trim_fastq_1():
    try:
        os.makedirs('test_out')
    except OSError:
        pass
    xout= oxbs_qc_func.trim_fastq('test_data/mjb042_oxBS_R1.txt.gz', 20, 'test_out/mjb042_oxBS_R1.fq')
    assert xout
    assert os.path.exists('test_out/mjb042_oxBS_R1.fq.gz') ## Note .gz included!

    ## Same number of lines in input and output
    tfq= gzip.open('test_out/mjb042_oxBS_R1.fq.gz').readlines()
    fq= gzip.open('test_data/mjb042_oxBS_R1.txt.gz').readlines()
    assert len(tfq) == len(fq)

    ## Lenght of read AND quality == 20
    read_quality_qlen= [len(tfq[x].strip()) for x in range(1, 100, 2)]
    assert set(read_quality_qlen) == set([20])

def test_task_bismark_exceptions():
    passed= False
    try:
        oxbs_qc_func.task_bismark_aln(inFastq= ['NA'], ref= 'control_reference/bsseq_synthetic4', outSam= 'mjb042_oxBS_R1.sam', opts= '') ## No input fastq
    except oxbs_qc_func.BismarkException:
        passed= True
    assert passed

    passed= False
    try:
        oxbs_qc_func.task_bismark_aln(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz'], ref= 'control_reference/bsseq_synthetic4', outSam= 'mjb042_oxBS_R1', opts= '') ## No exts on sam
    except oxbs_qc_func.BismarkException:
        passed= True
    assert passed

    passed= False
    try:
        oxbs_qc_func.task_bismark_aln(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz'], ref= 'control_reference/bsseq_synthetic4/', outSam= 'mjb042_oxBS_R1.sam', opts= '') ## No ref
    except oxbs_qc_func.BismarkException:
        passed= True
    assert passed

def test_task_bismark_aln_1():
    x= oxbs_qc_func.task_bismark_aln(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz'], ref= 'control_reference/bsseq_synthetic4', tmpdir= 'test_out', outSam= 'mjb042_oxBS_R1.sam', opts= '-o test_out')
    print(x)
    assert x['sam'] == 'test_out/mjb042_oxBS_R1.sam' ## Outut sam is what you expect
    assert os.path.exists(x['sam']) ## Sam exists
    assert len(x['all_out']) > 1 ## There are more files in output dir
    
def test_task_bismark_aln_2():
    x= oxbs_qc_func.task_bismark_aln(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz', 'test_data/mjb042_oxBS_R2.txt.gz'], ref= 'control_reference/bsseq_synthetic4', tmpdir= 'test_out', outSam= 'mjb042.sam', opts= '-o test_out')
    print(x)
    assert x['sam'] == 'test_out/mjb042.sam' ## Outut sam is what you expect
    assert os.path.exists(x['sam']) ## Sam exists
    assert len(x['all_out']) > 1 ## There are more files in output dir

def test_cleanSamReadNames():
    try:
        os.mkdir('test_out_cleanSamReadNames')
    except OSError:
        pass    
    shutil.copy('test_data/cc035.sam', 'test_out_cleanSamReadNames/cc035.sam')
    orifile= open('test_out_cleanSamReadNames/cc035.sam').readlines()
    oriStat= os.path.getmtime('test_out_cleanSamReadNames/cc035.sam')
    x= oxbs_qc_func.cleanSamReadNames('test_out_cleanSamReadNames/cc035.sam')
    newfile= open('test_out_cleanSamReadNames/cc035.sam').readlines()
    newStat= os.path.getmtime('test_out_cleanSamReadNames/cc035.sam')
    assert x
    assert len(orifile) == len(newfile)
    assert newStat != oriStat
    
def test_task_sam2bam():
    try:
        os.mkdir('test_out_sam2bam')
    except OSError:
        pass    
    shutil.copy('test_data/cc035.sam', 'test_out_sam2bam/cc035.sam')
    p= oxbs_qc_func.task_sam2bam('test_out_sam2bam/cc035.sam', rmSam= True, path= '')
    assert os.path.exists('test_out_sam2bam/cc035.bam')
    assert os.path.exists('test_out_sam2bam/cc035.bam.bai')
    assert p['bam'] == 'test_out_sam2bam/cc035.bam'
    shutil.rmtree('test_out_sam2bam')