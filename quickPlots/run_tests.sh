#!/usr/bin/env bash

## Move to repository dir containing plot and test_data/. The execute ./run_tests.sh

mkdir -p test_out/
rm test_out/test*

# echo -e "\n---- MAIN HELP ----\n"
./plot -h
./plot boxplot -h
./plot xyplot -h

# =============
 echo -e "\n---- HISTOGRAM ----\n"
 ./plot histogram -h
 ./plot histogram -i test_data/aln.sam -x 9 
 ./plot histogram -i test_data/aln.sam -x 9 -xwin 1
 ./plot histogram -i test_data/aln.sam -x 9 -b 10
 ./plot histogram -i test_data/aln.sam -x 9 -w 20 -h 10

 echo -e "\n---- BOXPLOT ----\n"
 ./plot boxplot -h
 ./plot boxplot -i test_data/aln.sam -x 2 -y 9
 ./plot boxplot -i test_data/aln.sam -x 2 -y 9 -bty varwidth
 ./plot boxplot -i test_data/aln.sam -x 2 -y 9 -bty notch

 echo -e "\n---- XYPLOT ----\n"
 ./plot xyplot -h
 ./plot xyplot -i test_data/aln.sam -x 4 -y 9 -s 2
 ./plot xyplot -i test_data/aln.sam -x 2 -y 5 -log xy
 ./plot xyplot -i test_data/cloud.txt -x 1 -y 2 -t smooth -d ','

echo -e "\n---- BARPLOT ----\n"
./plot barplot -h
./plot barplot -i test_data/bar.txt -x 2 -y 1
./plot barplot -i test_data/bar.txt -x 2 -y 1 -log y
./plot barplot -i test_data/bar.txt -x 2 -y 1 -ywin 1

echo -e "\n---- HEADER SELECTION ----\n"
# These two plots should be the same
./plot xyplot -i test_data/headers.txt -x 4 -y 2 -s 2 -d ','
./plot xyplot -i test_data/headers.txt -x x4 -y x2 -s 2 -d ','

## These should fail with col not found:
./plot xyplot -i test_data/headers.txt -x x4 -y na -s 2 -d ','
./plot xyplot -i test_data/headers.txt -x x4 -y 2 -s 2 -d ','

echo -e "\n---- CAN HANDLE MULTIPLE SEPARATORS ----"
awk '{print $1 "\t\t  \t\t" $2}' test_data/bar.txt | ./plot barplot -i - -x 2 -y 1 -d '\s+'
