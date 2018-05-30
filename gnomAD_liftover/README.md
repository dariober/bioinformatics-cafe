Liftover gnomAD to hg38
========================

[gnomAD](http://gnomad.broadinstitute.org/) files are based on hg19 with chromosome names in ENSEMBL format *i.e.*,
without `chr` prefix.

This script downloads, renames chromosomes and lifts over to hg38. In addition,
it removes all INFO tags other than AF. The final input is useful for *e.g.*
[gatk/getPileupSummaries](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.2.0/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.php).

See `liftover.smk` for exact commands and sources.

Note that at the time of this writing (30/05/2018), liftover files are available from gnomAD at:

```
gsutil ls gs://gnomad-public/release/2.0.2/vcf/genomes/liftover_grch38/
```

However, they are in ENSEMBL chromsome format and complete of all information (=they are big).

Usage
-----

```
snakemake --jobs 4 -p -s liftover.smk -C picard=/path/to/picard.jar ref=/path/to/hg38.genome.fa
```

To execute on cluster using `slurm`:

```
snakemake --jobs 100 -p -s liftover.smk -C picard=~/local/bin/picard.jar ref=/scratch/dberaldi/ref/GRCh38/GRCh38.fa \
    --cluster "sbatch --cpus-per-task={params.cpus} --mem={params.mem}"
```
