#!/bin/bash

# USAGE
# cd /.../bioinformatics-cafe/trunk/checkTerminatorBlock/test/
# bash run_tests.sh

echo -e "\nCAN SHOW HELP"
bash ../checkTerminatorBlock.sh -h

echo -e "\nCAN CHECK ONE FILES"
bash ../checkTerminatorBlock.sh hasblock.bam

echo -e "\nCAN GLOB FILES"
bash ../checkTerminatorBlock.sh *.bam

echo -e "\nMUST FAIL ON NON EXISTING FILE"
bash ../checkTerminatorBlock.sh not_exsisting.bam

echo -e "\nCAN PRINT EXIT CODE FROM LAST picard"
bash ../checkTerminatorBlock.sh hasblock.bam defective.bam > /dev/null
echo "Expected code 100. Actual: $?"
