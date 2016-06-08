#!/bin/bash

## Test against bedtools 24/05/2016 Version: v2.25.0

## Simple test
     fastaFromBed -fi test_data/ref.fa -bed test_data/a.bed -fo exp.tmp.fa
./fastaFromBed.py -fi test_data/ref.fa -bed test_data/a.bed -fo obs.tmp.fa
diff exp.tmp.fa obs.tmp.fa
rm exp.tmp.fa obs.tmp.fa

## Use name
     fastaFromBed -name -fi test_data/ref.fa -bed test_data/a.bed -fo exp.tmp.fa
./fastaFromBed.py -name -fi test_data/ref.fa -bed test_data/a.bed -fo obs.tmp.fa
diff exp.tmp.fa obs.tmp.fa
rm exp.tmp.fa obs.tmp.fa

## Use TAB
     fastaFromBed -tab -fi test_data/ref.fa -bed test_data/a.bed -fo exp.tmp.fa
./fastaFromBed.py -tab -fi test_data/ref.fa -bed test_data/a.bed -fo obs.tmp.fa
diff exp.tmp.fa obs.tmp.fa

## Use STRAND
     fastaFromBed -s -fi test_data/ref.fa -bed test_data/a.bed -fo exp.tmp.fa
./fastaFromBed.py -s -fi test_data/ref.fa -bed test_data/a.bed -fo obs.tmp.fa
diff exp.tmp.fa obs.tmp.fa

## All options combined
     fastaFromBed -name -tab -s -fi test_data/ref.fa -bed test_data/a.bed -fo exp.tmp.fa
./fastaFromBed.py -name -tab -s -fi test_data/ref.fa -bed test_data/a.bed -fo obs.tmp.fa
diff exp.tmp.fa obs.tmp.fa

## Zero length features: MUST FAIL
echo -e "chr1\t0\t0" |      fastaFromBed -fi test_data/ref.fa -bed - -fo exp.tmp.fa
echo -e "chr1\t0\t0" | ./fastaFromBed.py -fi test_data/ref.fa -bed - -fo obs.tmp.fa