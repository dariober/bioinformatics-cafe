#!/usr/bin/env bash

# USAGE
# ./test.sh

set -e 

source ../bashTestFunctions/bashTestFunctions.sh

mkdir -p tmp

pprint 'Split in equal chunks'

    rm -f tmp/*
    cat dat.txt | ./split.py -p tmp/x. -l 2
    n=`ls -1 tmp/* | wc -l`
    assertEquals 3 $n

    pprint 'Each file has n lines'
    n=`zcat tmp/x.0000.gz | wc -l`
    assertEquals 2 $n

    n=`zcat tmp/x.0001.gz | wc -l`
    assertEquals 2 $n

    n=`zcat tmp/x.0002.gz | wc -l`
    assertEquals 2 $n

pprint 'Reconstruct original input'

    exp=`md5sum dat.txt | cut -f 1 -d' '`
    obs=`zcat tmp/* | md5sum | cut -f 1 -d' '`
    assertEquals $exp $obs


pprint 'Last chunk is smaller'

    rm -f tmp/*
    cat dat.txt | ./split.py -p tmp/x. -l 4
    n=`ls -1 tmp/* | wc -l`
    assertEquals 2 $n

pprint 'Each file has n lines'

    n=`zcat tmp/x.0000.gz | wc -l`
    assertEquals 4 $n

    n=`zcat tmp/x.0001.gz | wc -l`
    assertEquals 2 $n

pprint 'Reconstruct original input'

    exp=`md5sum dat.txt | cut -f 1 -d' '`
    obs=`zcat tmp/* | md5sum | cut -f 1 -d' '`
    assertEquals $exp $obs

pprint 'nlines larger than file'

    rm -f tmp/*
    cat dat.txt | ./split.py -p tmp/x. -l 10
    n=`ls -1 tmp/* | wc -l` # Only one file created
    assertEquals 1 $n

pprint 'Each file has n lines'

    n=`zcat tmp/x.0000.gz | wc -l`
    assertEquals 6 $n

pprint 'Reconstruct original input'

    exp=`md5sum dat.txt | cut -f 1 -d' '`
    obs=`zcat tmp/* | md5sum | cut -f 1 -d' '`
    assertEquals $exp $obs

pprint 'Test empty input file'

    rm -f tmp/*
    cat empty.txt | ./split.py -p tmp/x.
    n=`ls -1 tmp/* | wc -l`
    assertEquals 1 $n
    n=`zcat tmp/x.0000.gz | wc -l`
    assertEquals 0 $n

rm -rf tmp
pprint 'ALL DONE'
