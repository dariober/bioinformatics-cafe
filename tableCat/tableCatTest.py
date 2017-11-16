#!/usr/bin/env py.test

"""
MEMO: Test functions MUST have "test_" as name prefix.
"""

import tableCat as tc
import pytest
import subprocess as sp

def test_globToList_canGlobFilesFromList():
    inglob= ['test/*.txt']
    globbed= tc.globToList(inglob)
    print globbed
    assert ['test/a.txt', 'test/b.txt', 'test/c.txt'] == globbed

    inglob= ['test/*.txt', 'test/*.txt', 'nonsense.txt']
    globbed= tc.globToList(inglob)
    assert ['test/a.txt', 'test/b.txt', 'test/c.txt'] == globbed

def test_xopen_canOpenFlatOrGzFile():
    fin= tc.xopen('test/a.txt')
    fin= tc.xopen('test/a.txt.gz')
    with pytest.raises(Exception):
        tc.xopen('nonsense')
    fin.close()

def test_parseID_canParseFilename():
    id= tc.parseID('test/a.txt', keepdir= False)
    assert 'a.txt' == id
    id= tc.parseID('test/a.txt', keepdir= True)
    assert 'test/a.txt' == id
    id= tc.parseID('test/a.txt', keepdir= False, regex= '\.txt$')
    assert 'a' == id

def test_runScript():
    cmd= './tableCat.py -h'
    p= sp.Popen(cmd, shell= True, stdout=sp.PIPE, stderr=sp.PIPE)
    print(p.stderr.read())
    p.communicate()
    assert p.returncode == 0
    
    cmd= './tableCat.py -i test/a.txt test/b.txt'
    p= sp.Popen(cmd, shell= True, stdout=sp.PIPE, stderr=sp.PIPE)
    print(p.stderr.read())
    p.communicate()
    assert p.returncode == 0

    cmd= './tableCat.py -i test/a.txt test/b.txt'
    p= sp.Popen(cmd, shell= True, stdout=sp.PIPE, stderr=sp.PIPE)
    out= p.stdout.read().split('\n')
    p.communicate()
    assert p.returncode == 0
    assert 'hdr1\thdr2\thdr3\ta.txt' == out[0]
    assert 'a\t1\t10\ta.txt' == out[1]

def test_runScriptReadFromStdin():
    cmd= 'ls test/a.txt test/b.txt | ./tableCat.py -i -'
    p= sp.Popen(cmd, shell= True, stdout=sp.PIPE, stderr=sp.PIPE)
    out= p.stdout.read().split('\n')
    p.communicate()
    assert p.returncode == 0
    assert 'hdr1\thdr2\thdr3\ta.txt' == out[0]
    assert 'a\t1\t10\ta.txt' == out[1]
