#!/bin/bash

python localEnrichmentBed.py -t test_data/fk135_BSF1_J.NPD_peaks.narrowPeak -b test_data/fk135_BSF1_J.tryp.bam -bl test_data/unmappable.150.bedgraph.gz > test.leb

exp=`wc -l test_data/fk135_BSF1_J.NPD_peaks.narrowPeak | awk '{print $1}'`
obs=`wc -l test.leb | awk '{print $1 -1}'` ## -1 to remove header

if [ "$exp" != "$obs" ]
then
    echo "WRONG: expected $exp got $obs"
else 
    echo "OK"
fi

##
exp=`echo 'Tb927_01_v5.1 34850 40660 8 66 232 39460 5810 236.68109812 4.57736776309' | tr ' ' '\t'`
obs=`awk 'NR == 9 {print $0}' test.leb`

if [ "$exp" != "$obs" ]
then
    echo -e "WRONG: expected\n $exp\n got\n $obs"
else 
    echo "OK"
fi

rm test.leb