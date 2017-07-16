Find a sequence motif in introns
================================

* Collect introns coordinates

These are the features in the GFF file:

```
cut -f3 data/PlasmoDB-30_PbergheiANKA.gff | grep -v '#'| sort | uniq -c
13799 CDS
13936 exon
5246 gene
5076 mRNA
  11 ncRNA
  53 rRNA
   5 snRNA
  31 snoRNA
  78 tRNA
   2 three_prime_UTR
```

Intron is every interval inside a *mRNA* or *gene* that doesn't overlap with any of the other features (*exons* and *CDS* mostly)

* Get coordinates of gene and mRNA

```
awk -v FS='\t' '$3 == "mRNA" || $3 == "gene"' data/PlasmoDB-30_PbergheiANKA.gff | sortBed | mergeBed > genes.bed
```

* Get coordinates of intra-gene features:

```
awk -v FS='\t' '$3 != "mRNA" && $3 != "gene"' data/PlasmoDB-30_PbergheiANKA.gff | sortBed | mergeBed > intragenes.bed
```

* Get introns:

```
subtractBed -a genes.bed -b intragenes.bed > introns.bed
```

* Collect coordinates of the motif

[fastaRegexFinder.py](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder)

`awk` strips the characters after the first blank space in the chromosome name.

```
fastaRegexFinder.py -f data/PlasmoDB-30_PbergheiANKA_Genome.fasta -r 'GTGTAC[GT]' \
| awk -v FS='\t' -v OFS='\t' '{gsub(" .*", "", $1); print $0}' \
| intersectBed -a - -b introns.bed -u > intron_motifs.bed
```

* Clean up

```
rm genes.bed intragenes.bed introns.bed 
```
