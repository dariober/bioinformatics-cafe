#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

echo "CAN SHOW HELP"
python ../sortBedAsBam.py -h 

echo "CAN SORT BY CHROM"
obs=`python ../sortBedAsBam.py -i unsorted.bed -b in.bam | cut -f4`
exp="f
e
g
i
h
a
b
d
c
z
zz
w"

if [[ $obs ==  $exp ]]
then
    printf "${GREEN}PASSED${NC}\n"
else
    printf "${RED}FAILED${NC}\n"
    echo -e "Expected:\n$exp"
    echo -e "Observed:\n$obs"
fi

echo "CAN SORT BY CHROM AND BY POSITION"
obs=`python ../sortBedAsBam.py -i unsorted.bed -b in.bam --sort | cut -f4`
exp="e
f
g
h
i
a
b
c
d
zz
z
w"

if [[ $obs ==  $exp ]]
then
    printf "${GREEN}PASSED${NC}\n"
else
    printf "${RED}FAILED${NC}\n"
    echo -e "Expected:\n$exp"
    echo -e "Observed:\n$obs"
fi