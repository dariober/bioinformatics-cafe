## Download fastq files from BaseSpace programmatically

For usage and examples see `BaseSpaceDownloader.R -h`.

## TODO

* ~~Check downloaded file size match size reported by BaseSpace~~ Done
* ~~If a file with the same name and size exists in the output dir, do not re-dowload it.~~ Done. However, it appears sometimes (?) the index files are downloaded first so they overwrite the existing files so this option gets useless! 
* Add a `--force` option to download even if file exists.
