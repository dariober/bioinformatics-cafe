# Extract methylation calls from bam file

Extract methylation calls from BAM file. Any bam file, sorted and indexed, produced by a bisulphite-aware aligner such as [bwameth](https://github.com/brentp/bwa-meth) should be suitable.

Help (*might be outdated*)

```
bam2methylation.py -h
usage: bam2methylation.py [-h] --input INPUT --ref REF [--l L] [--A]
                          [--samargs SAMARGS] [--region REGION [REGION ...]]
                          [--mismatch] [--minq MINQ] [--tmpdir TMPDIR]
                          [--keeptmp] [--quiet] [--version]

DESCRIPTION
    Extract methylation calls from BAM file.
    
OUTPUT:
   bedGraph with columns:
    <chrom>  <pos-1>  <pos>  <pct meth'd>  <cnt methylated>  <tot count>  <strand> [cnt mismatch]

SEE ALSO
https://github.com/dariober/bioinformatics-cafe/tree/master/bam2methylation

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Input bam file.
                                           
  --ref REF, -r REF     Reference fasta file.
                                           
  --l L, -l L           Bedfile with intervals where pileup should be generated.
                        Passed to `samtools mpileup -l <>`
                                           
  --A, -A               Passed to mpileup: Count anomalous read pairs. Default
                        is to exclude them.
                                           
  --samargs SAMARGS, -s SAMARGS
                        String of optional arguments passed to `samtools view` to filter
                        reads. Put this string in quotes leaving a space after opening quote. See also --region. 
                        E.g. -s ' -q 15 -F 256'
                                           
  --region REGION [REGION ...], -R REGION [REGION ...]
                        Region passed to `samtools view` to extract reads from.
                        E.g. -R 'chr1:0-10000'
                                           
  --mismatch, -mm       Insert a column of mismatches after the tot methylation count.
                                           
  --minq MINQ, -mq MINQ
                        Minimum base quality required to consider the base a
                        methylation call or a mismatch.
                                           
  --tmpdir TMPDIR       Temp directory. If not given python will find one.
                                           
  --keeptmp             Keep tmp dir. Use for debugging.
                                           
  --quiet               Do not print log messages
                                            
  --version             show program's version number and exit
```

## Usage

* Minimal example

```
bam2methylation.py -i aln.bam -r ref.fa

chr10   101082  101083  25.00   1       4       -
chr10   101086  101087  25.00   1       4       -
chr10   101119  101120  0.00    0       3       -
chr10   101325  101326  0.00    0       3       +
chr10   101326  101327  0.00    0       1       -
chr10   102013  102014  50.00   1       2       +
...
```

* Return only target positions and filter out reads with mapq < 10, exclude reads with flag 256, and ignore base quality < 15

```
bam2methylation.py -i aln.bam -r ref.fa -l targets.bed -s ' -q 15 -F 256' -mq 15
```

* Get only the Cs in a certain context (e.g. CpGs)

You need to pass a bed file where these positions are. For example, first
extract all the CpGs in the genome with
[fastaRegexFinder.py](https://raw.githubusercontent.com/dariober
/bioinformatics-cafe/master/fastaRegexFinder.py):

```
fastaRegexFinder.py -f mmu.fa -r CG --noreverse > mm9.allcpg.bed
```

Then pass the bed file of regions of interest to `-l <reg.bed>` option.

* Dumb parallelisation

A simple way to speed up the process: Process each chromosome separately using the -R option:

```
for chr in chr1 chr2 chr3 ...
do
    echo "bam2methylation.py -i aln.bam -r ref.fa -R $chr > aln.${chr}.tmp.bedgraph" > aln.${chr}.tmp.sh
done

## Execute in parallel using e.g. xargs, parallel, bsub or else:
ls aln.*.tmp.sh | xargs -P 12 -n 1 bash && rm aln.*.tmp.sh

## Once done, concatenate individual chromosomes:
cat aln.*.tmp.bedgraph > aln.bedgraph && rm aln.*.tmp.bedgraph
```

## Output format

Output goes to `stdout` as tab delimited columns in extended bedGraph format with columns:

```
<chrom>  <pos-1>  <pos>  <pct meth'd>  <cnt methylated>  <tot count>  <strand> [cnt mismatch]
```

**Memo**: bedGraph is 0-based, so if the first base of chr1 is C it will have position: `chrom 0 1 ... +`. Output is sorted by chromosome and position.

With `--mismatch` option an additional column of mismatch count printed as last column. This is the number of A or G when the reference is C.

See also [addContextToBed.sh](https://github.com/dariober/bioinformatics-cafe/tree/master/addContextToBed) to add sequence context to each position.

## Installation and input requirements

Download, make executable and, optionally, copy to directory on PATH (*e.g.* `~/bin/`):

```
wget https://raw.githubusercontent.com/dariober/bioinformatics-cafe/master/bam2methylation/bam2methylation.py
chmod a+x bam2methylation.py
mv bam2methylation.py /dir/on/path/bin/
```

Requirments:

* Python 2.7 (python 2.6 or 3.x should also work)
* [samtools](http://www.htslib.org/) version 1.1+ on `PATH`
* Unix `sort` and `awk`

Required input is a bam file sorted and indexed and the corresponding reference sequence in fasta format.
