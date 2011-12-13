## Align CAGE tags hg95ctrls from human macrophages against human genome hg19. 
## Return one match for read with multiple hits (-r 1 option).

/exports/work/vet_roslin_nextgen/SOAP/soap2.19release/soap \
  -a /exports/work/vet_roslin_nextgen/dario/fastq/hg95ctrls.fa \
  -D /exports/work/vet_roslin_nextgen/dario/soap/indexfiles/Hsapiens/hg19.fa.index \
  -o /exports/work/vet_roslin_nextgen/dario/soap/output/aligned/hg95ctrls/20091019_hg95ctrls_vs_Human_hg19_r1.map \
  -r 1

