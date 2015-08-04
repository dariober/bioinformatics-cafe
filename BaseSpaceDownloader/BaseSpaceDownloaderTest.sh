#!/bin/bash

# This project ID referes to Mike Booth's run:
# Contains run https://basespace.illumina.com/run/2365367/MspI-enrich-test
# If this project or run is changed the tests below will fail!!
TEST_PROJ=1451454 

bsd=$1
if [[ $bsd == '' ]]
then
echo ''
echo 'USAGE'
echo 'BaseSpaceDownloaderTest.sh /path/to/BaseSpaceDownloader.R'
echo 'All test should return PASSED'
exit
fi
# ==============================================================================

echo -e '\nCAN SHOW HELP'
$bsd -h

echo -e '\nCAN echo SAMPLES'
# ====================
$bsd -p $TEST_PROJ -e > out.txt
result=`diff tests/expect-echo.txt out.txt`
if [[ $result == '' ]]
then
    echo 'PASSED'
else
    echo '**** FAILED ****'
    echo $results
fi
rm out.txt

echo -e '\nCAN GREP FILENAME'
# ======================
$bsd -p $TEST_PROJ -r 'MspI-[AC]' -e > out.txt
result=`diff tests/expect-grep-filename.txt out.txt`
if [[ $result == '' ]]
then
    echo 'PASSED'
else
    echo '**** FAILED ****'
fi
rm out.txt

echo -e '\nCAN GREP SAMPLE'
# ======================
$bsd -p $TEST_PROJ -s 'MspI-[DF]' -e > out.txt
result=`diff tests/expect-grep-sample.txt out.txt`
if [[ $result == '' ]]
then
    echo 'PASSED'
else
    echo '**** FAILED ****'
fi
rm out.txt

echo -e '\nCAN DOWNLOAD SAMPLES (can take ages)'
# ========================
$bsd -p $TEST_PROJ -r 'MspI-[DF].*_L001_R1_001.fastq.gz'

test -e Data/Intensities/BaseCalls/MspI-D_S4_L001_R1_001.fastq.gz
if [[ $? == 0 ]]
then
    echo 'PASSED'
else
    echo '**** FAILED ****'
fi

test -e Data/Intensities/BaseCalls/MspI-F_S6_L001_R1_001.fastq.gz 
if [[ $? == 0 ]]
then
    echo 'PASSED'
else
    echo '**** FAILED ****'
fi
rm -i -r Data