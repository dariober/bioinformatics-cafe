#!/use/bin/env bash

n=`./catInterleaveFastq.sh -1 data/r1.fq -2 data/r2.fq | wc -l`
if [[ $n != 24 ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

n=`./catInterleaveFastq.sh -1 data/r1.fq data/r1.fq -2 data/r2.fq data/r2.fq | wc -l`
if [[ $n != 48 ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

x=`./catInterleaveFastq.sh -1 data/r1.fq data/r1.fq -2 data/r2.fq data/r2.fq | head -n 1`
if [[ $x !=  '@r1.1' ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

x=`./catInterleaveFastq.sh -1 data/r1.fq data/r1.fq -2 data/r2.fq data/r2.fq | tail -n+5 | head -n 1`
if [[ $x !=  '@r1.2' ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

# Name of last read
x=`./catInterleaveFastq.sh -1 data/r1.fq data/r1.fq -2 data/r2.fq data/r2.fq | tail -4 | head -n 1`
if [[ $x !=  '@r3.2' ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

# Order of -1 and -2 doesn't matter
x=`./catInterleaveFastq.sh -2 data/r2.fq data/r2.fq -1 data/r1.fq data/r1.fq | head -n 1`
if [[ $x !=  '@r1.1' ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

n=`./catInterleaveFastq.sh -1 data/r1.fq.gz data/r1.fq.gz -2 data/r2.fq.gz data/r2.fq.gz | wc -l`
if [[ $n != 48 ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

# Unequal number of files
./catInterleaveFastq.sh -2 data/r1.fq data/r1.fq -1 data/r2.fq > /dev/null
if [[ $? == 0 ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

# Positional argument 2 not given
./catInterleaveFastq.sh -1 data/r1.fq data/r1.fq > /dev/null
if [[ $? == 0 ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

# Print version & exit
version=`./catInterleaveFastq.sh -v`
if [[ $? != 0 || -z $version ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

echo "All OK"
