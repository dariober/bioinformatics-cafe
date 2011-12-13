#!/bin/bash

## Submits all the following scripts to eddie in such way that one script
## doesn't start until the previous one has completed.

# Move to working dir
cd /exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment

cp 20101218_merge_cage_sam_iterated.py merge_cage_sam_iterated.py
qsub -P vet_roslin -pe memory 1 -l h_rt=00:29:00 -hold_jid bwa_27.sh merge_cage_sam_iterated.py

cp 20101218_rDNA_filter.py rDNA_filter.py
qsub -P vet_roslin -pe memory 1 -l h_rt=00:29:00 -hold_jid merge_cage_sam_iterated.py rDNA_filter.py

cp 20101218_sam_single_mappers.py sam_single_mappers.py
qsub -P vet_roslin -pe memory 1 -l h_rt=00:29:00 -hold_jid rDNA_filter.py sam_single_mappers.py

cp 20101218_sam_to_bam.sh sam_to_bam.sh
qsub -P vet_roslin -pe memory 1 -l h_rt=00:29:00 -hold_jid sam_single_mappers.py sam_to_bam.sh

cp 20101218_pileup.sh pileup.sh
qsub -P vet_roslin -pe memory 1 -l h_rt=00:29:00 -hold_jid sam_to_bam.sh pileup.sh

## Rememebr to change the input for this manually!!!
cp 20101218_parse_pileup.py parse_pileup.py
qsub -P vet_roslin -pe memory 1 -l h_rt=00:29:00 -hold_jid pileup.sh parse_pileup.py