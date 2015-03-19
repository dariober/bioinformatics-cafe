#! /usr/bin/env Rscript

library('RUnit')

# Usage:
# ------------------------------------------------------------------------------
# Add tests to file tests/<foo>Test.R
# Then execute
# ./run_tests.R
# run_tests.R itself shouldn't need edits.
# ------------------------------------------------------------------------------
# This is an awful workaround to load only the functions w/o executing the argparse
# bit.
# Similar to python if __name__ == "__main__"
sc<- pipe("awk '{if($0 ~ /END_OF_FUNCTIONS/) {exit}else{print $0}}' BaseSpaceDownloader.R")
source(sc)
close(sc)

# source('BaseSpaceDownloader.R') 
 
test.suite <- defineTestSuite("TestSuite.BaseSpaceDownloader",
                              dirs = file.path("tests"),
                              testFileRegexp = 'Test\\.R$')
 
test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)