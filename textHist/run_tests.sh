#!/usr/bin/env bash

echo "CAN PRINT HELP"
./textHist.R -h

echo "CAN PLOT HISTOGRAM DEFAULT"
./textHist.R -i test_data/data.txt -c 1 

echo "CAN PLOT COL WITH MISSING DATA"
./textHist.R -i test_data/data.txt -c 2 -V

echo "CAN SET WIDTH"
./textHist.R -i test_data/data.txt -c 1 -w 80
