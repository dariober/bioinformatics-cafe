#!/bin/bash

# =========================  SET UP ===========================================

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

PASSED="${GREEN}PASSED${NC}\n"
FAILED="${RED}FAILED${NC}\n"

tmp='tmpRunTests'
mkdir $tmp

# =========================  TESTS ===========================================

echo "CAN CENTER INTERVALS"
../makeTargetBed.sh -g test_data/genome.txt -i test_data/target.bed --center --slop 0 --binSize 1 > $tmp/obs_centerInt.bed
result=`diff test_data/exp_centerInt.bed $tmp/obs_centerInt.bed`

if [[ $result == '' ]]
then
    printf $PASSED
else
    printf $FAILED
fi

echo "CAN MAKE WINDOWS"
../makeTargetBed.sh -g test_data/genome.txt -i test_data/target.bed -c -s 3 -b 1 > $tmp/obs_windows.bed
result=`diff test_data/exp_windows.bed $tmp/obs_windows.bed`

if [[ $result == '' ]]
then
    printf $PASSED
else
    printf $FAILED
    echo "Expected"
    cat test_data/exp_windows.bed
    echo "Observed"
    cat $tmp/obs_windows.bed
fi

echo "CAN READ STDIN"
cat test_data/target.bed | ../makeTargetBed.sh -g test_data/genome.txt -i - -c -s 3 -b 1 2> $tmp/std.err

if [[ $? == 0 ]]
then
    printf $PASSED
else
    printf $FAILED
    cat $tmp/std.err
fi

# ------------------------------------------------------------------------------
# ../makeTargetBed.py -g test_data/genome.txt -i test_data/target.bed -c -s 3 -b 1 > test_data/targetWindows.bed
echo "CAN MAKE PROFILE"
../bedProfiler.py -R test_data/targetWindows.bed -s test_data/scores.bedGraph \
| awk '{if(NR > 1){printf "%s\t%s\t%.1f\t%s\t%s\n", $1, $2, $3, $4, $5} else {print $0}}' > $tmp/obs_profile.tab
result=`diff test_data/exp_profile.tab $tmp/obs_profile.tab`

if [[ $result == '' ]]
then
    printf $PASSED
else
    printf $FAILED
    echo "Expected"
    cat test_data/exp_profile.tab
    echo "Observed"
    cat $tmp/obs_profile.tab
fi

# ========================== CLEAN UP ==========================================

rm -r $tmp/
