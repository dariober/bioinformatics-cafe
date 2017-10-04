Interleave fastq files
======================

This script interleaves one or more pairs of fastq files. 
It's useful when several pairs of fastq files from the same library needs to be processed by a downstream programs, typically
an aligner or an adapter trimmer that accepts interleaved input from `stdin`:

```
catInterleaveFastq.sh sampleA_R1.L001.fq.gz,sampleA_R1.L002.fq.gz \
                      sampleA_R2.L001.fq.gz,sampleA_R2.L002.fq.gz \
| bwa mem -p genome.fa - \
| ...
```

Download and run
----------------

To download only this script:

```
curl -O https://raw.githubusercontent.com/dariober/bioinformatics-cafe/master/catInterleaveFastq.sh
chmod 755 catInterleaveFastq.sh
./catInterleaveFastq.sh -h # Show help
```
