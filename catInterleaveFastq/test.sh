#!/use/bin/env bash

n=`./catInterleaveFastq.sh data/r1.fq data/r2.fq | wc -l`
if [[ $n != 24 ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

n=`./catInterleaveFastq.sh data/r1.fq,data/r1.fq data/r2.fq,data/r2.fq | wc -l`
if [[ $n != 48 ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

x=`./catInterleaveFastq.sh data/r1.fq,data/r1.fq data/r2.fq,data/r2.fq | head -n 1`
if [[ $x !=  '@r1.1' ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

x=`./catInterleaveFastq.sh data/r1.fq,data/r1.fq data/r2.fq,data/r2.fq | tail -n+5 | head -n 1`
if [[ $x !=  '@r1.2' ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

n=`./catInterleaveFastq.sh data/r1.fq.gz,data/r1.fq.gz data/r2.fq.gz,data/r2.fq.gz | wc -l`
if [[ $n != 48 ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

# Unequal number of files
./catInterleaveFastq.sh data/r1.fq,data/r1.fq data/r2.fq > /dev/null
if [[ $? == 0 ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

# Positional argument 2 not given
./catInterleaveFastq.sh data/r1.fq,data/r1.fq > /dev/null
if [[ $? == 0 ]]
then
    echo "WRONG at line $LINENO"
    exit 1
fi

echo "All OK"
