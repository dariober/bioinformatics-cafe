## Align CAGE tags hg95ctrls from human macrophages against pig genome. 
## Return one match for read with multiple hits (-r 1 option).

/exports/work/vet_roslin_nextgen/SOAP/soap2.19release/soap \
  -a /exports/work/vet_roslin_nextgen/dario/fastq/hg95ctrls.fa \
  -D /exports/work/vet_roslin_nextgen/dario/soap/indexfiles/Hsapiens/Human_ribosomal_DNA_complete_repeating_unit.fa.index \
  -o /exports/work/vet_roslin_nextgen/dario/soap/output/aligned/20091014_hg95ctrls_vs_Human_rDNA_r1.map \
  -r 1

