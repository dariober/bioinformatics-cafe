Quantification of read enrichment over flanking background
==========================================================

Compute the read enrichment in target intervals relative to local background.

Typical use case: A ChIP-Seq experiment on a sample returns a number of regions
of enrichment. We want to know how enriched these regions are in a *different*
sample. Note that enrichment is quantified relative to the local background
not relative to an input control.

Getting started
---------------

Download, make executable, show help:

```
wget 'https://github.com/dariober/bioinformatics-cafe/blob/master/localEnrichmentBed/localEnrichmentBed.py?raw=true' -O localEnrichmentBed.py

chmod a+x localEnrichmentBed.py

./localEnrichmentBed.py -h
```


See also on Biostars https://www.biostars.org/p/195689/