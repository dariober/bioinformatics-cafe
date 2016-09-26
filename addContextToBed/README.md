Description
===========

Add sequence context to intervals in a bed file

Usage
======

```
addContextToBed.sh <bedgraph> <fasta-ref> <strand-column-index> <up> <down>
```

Positional arguments

1. Bed file to annotate 
2. Reference fasta file. *.fai index must be present in the same dir.
3. Column index (1-based) in bed file with strand info. '+' for forward '-' for reverse
4. Add this many bases to from UPstream of the interval (i.e. to the left if +strand on forw)
5. Add this many bases to from DOWNstream of the C (i.e. to the right if -strand on forw)


Requirements
============

* [bedtools](http://bedtools.readthedocs.io/en/latest/) on path

* [fastaFromBed.py](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaFromBed) on path (`bedtools fastaFromBed` 
   will also work but it will be very slow for long files)

* Fasta file indexed with index named <seqname>.fai in the same dir as input fasta. You can use `samtools faidx` to create it.

  
Notes
=====

* Error handling is not very robust, if one of the pipes fails the script keeps going!

* If an extended interval is beyond a chromosome length, the reported sequence is truncated.
