#!/bin/bash

# Run test by running 
# bash run_tests.sh

bash genomeFileFromBam.sh test_data/test.bam > test_data/obs.genome

diff test_data/obs.genome test_data/exp.genome > test_data/diff.tmp

if [ ! -s test_data/diff.tmp ]; then
    echo "TEST OK"
    rm test_data/diff.tmp
    rm test_data/obs.genome
else
    cat test_data/diff.tmp
fi

