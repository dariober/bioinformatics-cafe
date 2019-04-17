#!/bin/bash

python localEnrichmentBed.py -t test_data/fk135_BSF1_J.NPD_peaks.narrowPeak -b test_data/fk135_BSF1_J.tryp.bam -bl test_data/unmappable.150.bedgraph.gz > test.leb

exp=`wc -l test_data/fk135_BSF1_J.NPD_peaks.narrowPeak | awk '{print $1}'`
obs=`wc -l test.leb | awk '{print $1}'` 

if [ "$exp" != "$obs" ]
then
    echo "WRONG: expected $exp got $obs"
    exit 1
else 
    echo "OK"
fi

##
exp=`echo 'Tb927_01_v5.1 34850 40660 fk135_BSF1_J.NPD_peak_8 29075 . 7.52465 2907.52271 2903.27173 2227 66 232 39460 5810 236.681 4.576' | tr ' ' '\t'`
obs=`awk 'NR == 8 {print $0}' test.leb`

if [ "$exp" != "$obs" ]
then
    echo -e "WRONG: expected\n $exp\n got\n $obs"
    exit 1
else 
    echo "OK"
fi

# Read stdin with minimal columns

cut -f 1-3 test_data/fk135_BSF1_J.NPD_peaks.narrowPeak \
| python localEnrichmentBed.py -t - -b test_data/fk135_BSF1_J.tryp.bam > test.leb

exp=`wc -l test_data/fk135_BSF1_J.NPD_peaks.narrowPeak | awk '{print $1}'`
obs=`wc -l test.leb | awk '{print $1}'` 

if [ "$exp" != "$obs" ]
then
    echo "WRONG: expected $exp got $obs"
    exit 1
else 
    echo "OK"
fi

# No records in input is ok

head -n 0 test_data/fk135_BSF1_J.NPD_peaks.narrowPeak \
| python localEnrichmentBed.py -t - -b test_data/fk135_BSF1_J.tryp.bam > test.leb

if [ "$?" != "0" ]
then
    echo "WRONG: expected exit code 0"
    exit 1
else 
    echo "OK"
fi

rm test.leb
