#!/bin/bash

PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/fastx_toolkit/fastx_toolkit_0.0.13/bin

fastx_quality_stats \
    -i /exports/work/vet_roslin_nextgen/dario/fastq/20100805_CAGE_AM/Beraldi_CAGE_50810_101202_EBRI093151_0060_s_8_sequence.txt \
    -o 20101206_CAGE_qc.txt