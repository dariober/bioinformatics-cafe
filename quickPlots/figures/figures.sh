# Script for figures to put in docs
## cd to ./plot script
# cd /Users/berald01/svn_git/bioinformatics-cafe/trunk/quickPlots/

# HISTOGRAM EXAMPLE
./plot histogram -i test_data/aln.sam -x 9 -o figures/hist.png -w 14 -h 10
./plot histogram -i test_data/aln.sam -x 9 -xwin 2 -o figures/hist-win.png -w 14 -h 10

./plot boxplot -i test_data/data.txt -x strain -y height -o figures/boxplot.png -w 12 -h 10
./plot boxplot -i test_data/data.txt -x strain -y height -ywin 2 -o figures/boxplot-ywin.png -w 12 -h 10
./plot boxplot -i test_data/data.txt -x strain -y height -ywin 2 -bty violin -o figures/violin-ywin.png -w 12 -h 10

./plot xyplot test_data/cloud.txt -i test_data/cloud.txt -x height -y weight -p 4 -d ',' -o figures/xyplot-1.png -w 12 -h 12
./plot xyplot test_data/cloud.txt -i test_data/cloud.txt -x height -y weight -p 4 -a 0.4 -d ',' -o figures/xyplot-2.png -w 12 -h 12
./plot xyplot test_data/cloud.txt -i test_data/cloud.txt -x height -y weight -t smooth -d ',' -o figures/smooth.png -w 12 -h 12