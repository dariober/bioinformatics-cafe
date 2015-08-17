#!/bin/bash

# USAGE
# cd /.../bioinformatics-cafe/trunk/checkTerminatorBlock/test/
# bash run_tests.sh

echo -e "\nCAN SHOW HELP"
../checkTerminatorBlock.sh -h

echo -e "\nCAN CHECK ONE FILES"
../checkTerminatorBlock.sh hasblock.bam

echo -e "\nCAN GLOB FILES"
../checkTerminatorBlock.sh *.bam

echo -e "\nCAN FAIL ON NON EXISTING FILE"
../checkTerminatorBlock.sh not_exsisting.bam