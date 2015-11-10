#!/bin/bash

# MEMO: If you change ref.fa update the index with `samtools faidx ref.fa`
# USAGE: Should return PASS
#  ./run_tests.sh

out=`../addContextToBed.sh in.bed ref.fa 7 2 3 | awk 'NR == 1 {print $NF}'`
if [[ $out == 'CTAC' ]]; then
    echo "PASS"
else
    echo "*** FAILED ***"
fi

out=`../addContextToBed.sh in.bed ref.fa 7 2 3 | awk 'NR == 2{print $NF}'`
if [[ $out == 'TACAATNN' ]]; then
    echo "PASS"
else
    echo "*** FAILED ***"
fi

out=`../addContextToBed.sh in.bed ref.fa 7 2 3 | awk 'NR == 3{print $NF}'`
if [[ $out == 'AACNNA' ]]; then
    echo "PASS"
else
    echo "*** FAILED ***"
fi

# Gzip files
out=`../addContextToBed.sh in.bed.gz ref.fa 7 2 3 | awk 'NR == 1 {print $NF}'`
if [[ $out == 'CTAC' ]]; then
  echo "PASS"
else 
  echo "*** FAILED ***"
fi
