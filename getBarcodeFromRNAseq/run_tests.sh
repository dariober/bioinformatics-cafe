#/bin/env bash

python getBarcodeFromRNAseq.py \
    -l 5 \
    -inR1 test_data/r1.fq \
    -inR2 test_data/r2.fq \
    -outR1 test_data/obs1.fq \
    -outR2 test_data/obs2.fq

n=`head -n 1 test_data/obs1.fq | grep -P '^@TTAAT:' | wc -l`
if [[ $n != 1 ]]
then
    echo "WRONG 1"
fi

n=`head -n 1 test_data/obs2.fq | grep -P '^@TTAAT:' | wc -l`
if [[ $n != 1 ]]
then
    echo "WRONG 2"
fi

rm test_data/obs1.fq test_data/obs2.fq

## Check we are not leaving behind any char
python getBarcodeFromRNAseq.py \
    -l 1 \
    -s '' \
    -inR1 test_data/r1.fq \
    -inR2 test_data/r2.fq \
    -outR1 test_data/obs1.fq \
    -outR2 test_data/obs2.fq

expectedWc=`cat test_data/r1.fq test_data/r2.fq | wc`
observedWc=`cat test_data/obs1.fq test_data/obs2.fq | wc`

if [[ "$expectedWc" != "$observedWc" ]]
then
    echo "WRONG 3"
fi

## test stdin stdout
python getBarcodeFromRNAseq.py \
    -inR1 <(cat test_data/r1.fq) \
    -inR2 <(cat test_data/r2.fq) | wc
