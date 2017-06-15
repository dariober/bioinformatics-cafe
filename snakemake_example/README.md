Example snakemake pipeline
==========================

```
snakemake --printshellcmds --dryrun -s scripts/cramToCounts.snake \
    --config cram=data/cram/*.cram \
             REF=data/static/genome.fa \
             GFF=data/static/genes.gff 
```

`config` arguments
------------------

* `cram`: Wild card (glob in fact) to capture the cram files to process.

* `REF`: Reference sequence for alignment.

* `GFF`: Reference gff file for alignment and count of reads in genes.
