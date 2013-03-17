#!/usr/bin/env py.test

import fastqc_md5
import shlex
import subprocess
import tempfile
import os
import shutil
from distutils import spawn

def test_getFastqcOutdir():
    fastqccmd= shlex.split('-t 8 --noextract -o /test/dir')    
    outd= fastqc_md5.getFastqcOutdir(fastqccmd)
    assert outd == '/test/dir'

def test_md5sum():
    tmpf= tempfile.NamedTemporaryFile()
    md5= fastqc_md5.md5sum(tmpf.name)
    assert md5 == 'd41d8cd98f00b204e9800998ecf8427e' ## This is md5 of NULL. possibly make better test

def test_add_md5_fastqc():
    tmpf= tempfile.NamedTemporaryFile(prefix= 'fastqc_data_', delete= True)
    shutil.copyfile('fastqc_data.txt', tmpf.name)
    fastqcList= fastqc_md5.add_md5_fastqc(tmpf.name, 'd41d8cd98f00b204e9800998ecf8427e')
    assert 'md5sum\td41d8cd98f00b204e9800998ecf8427e\n' in fastqcList

def test_false_if_fastqc_not_available():
    fastqc= fastqc_md5.fastqc_available('non/existant/path')
    assert fastqc is False

def test_fastqc_on_path():
    fastq_on_path= spawn.find_executable('fastqc')
    assert fastq_on_path is not None

def test_run():
    cmd= "./fastqc_md5.py -i db001.test.fq.gz --fastqc ' --noextract ' "
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stderr, stdout= p.communicate()
    done= True
    try:
        os.remove('db001.test.fq_fastqc.zip')
    except OSError:
        done= False
    assert done
    
    
    