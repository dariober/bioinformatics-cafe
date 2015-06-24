#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

echo "CROSS TAB TWO FILES"
../crossTabulate.py -n 1 -c 3 file01.txt file02.txt > test.out

result=`diff expected01.out test.out`
if [[ $result == '' ]]
then
printf "${GREEN}PASSED${NC}\n"
else
printf "${RED}FAILED${NC}\n"
fi
rm test.out

echo "FAIL IF NAMES ARE DIFFERENT (You should now see an error)"
../crossTabulate.py -n 1 -c 3 file01.txt file03.txt > /dev/null
if [[ 1 == $?  ]]
then
printf "${GREEN}PASSED${NC}\n"
else
printf "${RED}FAILED${NC}\n"
fi

echo "CAN GLOB FILES"
../crossTabulate.py -n 1 -c 3 file0[1-2].txt > test.out

result=`diff expected01.out test.out`
if [[ $result == '' ]]
then
printf "${GREEN}PASSED${NC}\n"
else
printf "${RED}FAILED${NC}\n"
fi
rm test.out

echo "CAN STRIP REGEX FROM FILE NAMES"
../crossTabulate.py -s '\.txt$' -n 1 -c 3 file0[1-2].txt > test.out

result=`diff expected02.out test.out`
if [[ $result == '' ]]
then
printf "${GREEN}PASSED${NC}\n"
else
printf "${RED}FAILED${NC}\n"
fi
rm test.out

echo "CAN PROCESS GZIP FILES"
../crossTabulate.py -n 1 -c 3 file01.txt.gz file02.txt.gz > test.out

result=`diff expected03.out test.out`
if [[ $result == '' ]]
then
printf "${GREEN}PASSED${NC}\n"
else
printf "${RED}FAILED${NC}\n"
fi
rm test.out

