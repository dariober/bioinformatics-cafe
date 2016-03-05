## Extract methylation calls from bam file

Extract methylation calls from BAM file. Only positions with non-zero count are printed.
If you want to get only the Cs in a certain context (e.g. CpGs) you need to pass
a bed file where these positions are. For example, first extract all the CpGs in the genome with
[fastaRegexFinder.py](https://raw.githubusercontent.com/dariober/bioinformatics-cafe/master/fastaRegexFinder.py):

```
fastaRegexFinder.py -f mmu.fa -r CG --noreverse > mm9.allcpg.bed
```

Then pass the bed file of regions of interest to `-l <reg.bed>` option.

## Usage and output

Minimal example

```
bam2methylation.py -i aln.bam -r ref.fa
```

Return only target positions and filter out reads with mapq < 10, exclude reads with flag 256, and ignore base quality < 15

```
bam2methylation.py -i aln.bam -r ref.fa -l targets.bed -s ' -q 15 -F 256' -mq 15
```

Output goes to `stdout` as tab delimited columns in extended bedGraph format with columns:

```
<chrom>  <pos-1>  <pos>  <pct meth'd>  <cnt methylated>  <tot count>  <strand> [cnt mismatch]
```

**Memo**: bedGraph is 0-based, so if the first base of chr1 is C it will have position: `chrom 0 1 ... +`. Output is sorted by chromosome and position.

With `--mismatch` option an additional column of mismatch count printed as last column. This is the number of A or G when the reference is C.

See also [addContextToBed.sh](https://github.com/dariober/bioinformatics-cafe/tree/master/addContextToBed) to add sequence context to 
each position.

## Installation and input requirements

Download, make executable and, optionally, copy to directory on PATH (**e.g.** `~/bin/`):

```
wget https://raw.githubusercontent.com/dariober/bioinformatics-cafe/master/bam2methylation/bam2methylation.py
chmod a+x bam2methylation.py
mv bam2methylation.py /dir/on/path/bin/
```

Requirments:

* Python 2.7 (python 2.6 or 3.x should also work)
* [samtools](http://www.htslib.org/) on `PATH`
* Unix `sort` and `awk`

Required input is a bam file sorted and indexed and the corresponding reference sequence in fasta format.
