#!/bin/bash

bash interleaveFastq.sh test_data/fq1.fq test_data/fq2.fq > test_data/obs.fq
xdiff=`diff test_data/expected.fq test_data/obs.fq`
if [[ "$xdiff" != "" ]]
then
    echo "WRONG 1"
    exit 1
fi

bash interleaveFastq.sh test_data/fq1.fq.gz test_data/fq2.fq.gz > test_data/obs.fq
xdiff=`diff test_data/expected.fq test_data/obs.fq`
if [[ "$xdiff" != "" ]]
then
    echo "WRONG 2"
    exit 1
fi

bash interleaveFastq.sh test_data/fq1.fq test_data/fq2.fq.gz > test_data/obs.fq
xdiff=`diff test_data/expected.fq test_data/obs.fq`
if [[ "$xdiff" != "" ]]
then
    echo "WRONG 3"
    exit 1
fi

bash interleaveFastq.sh test_data/fq1.fq.gz test_data/fq2.fq > test_data/obs.fq
xdiff=`diff test_data/expected.fq test_data/obs.fq`
if [[ "$xdiff" != "" ]]
then
    echo "WRONG 4"
    exit 1
fi

## Reading from stdin
bash interleaveFastq.sh <(cat test_data/fq1.fq) <(cat test_data/fq2.fq) > test_data/obs.fq
xdiff=`diff test_data/expected.fq test_data/obs.fq`
if [[ "$xdiff" != "" ]]
then
    echo "WRONG 5"
    exit 1
fi

bash interleaveFastq.sh test_data/fq1.fq.bz2 test_data/fq2.fq.bz2 > test_data/obs.fq
xdiff=`diff test_data/expected.fq test_data/obs.fq`
if [[ "$xdiff" != "" ]]
then
    echo "WRONG 6"
    exit 1
fi

rm test_data/obs.fq
echo "All tests passed"