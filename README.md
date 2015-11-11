# bioinformatics-cafe
Unsorted scripts for bioinformatics.

Formerly was on google.code as [bioinformatics-misc](https://bioinformatics-misc.googlecode.com/)

test commit

# TODO
`bam2methylation.py` 
* Replace individual options passed to `mpileup` (`-A, -l` etc) with an option that takes a string of args passed to `mpileup` as is. E.g. `--margs "-A -l regions.bed -d 100000"`.
* Remove option `--minq`. `mpileup` has `-Q` to do that same thing!
* Remove `--samargs` and calls to `samtools view`. As per `samtools v1.1`, `mpileup` options should replace the call top `view` completely.  For example, to include/exclude reads, instead of `samtools view -F/-f` use `samtools mpileup --ff/--rf`. 
