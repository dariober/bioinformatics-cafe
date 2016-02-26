## A command line utility for quick data exploration

Pourpose: Produce various types of plots directly from command line in order to get a quick overview of the input data. There should be minimal need to tweak input data and plot parameters. Plots are rendered using [R/ggplot](http://ggplot2.org/) but there's no need to know anything about R or ggplot.

### Use cases

* What is the distribution of insert size in a bam file? Insert size is in column 9:

```
samtools view aln.bam | plot histogram -i - -x 9
```

To avoid reading in the entire file let's take every *N*th read:

```
samtools view aln.bam | awk 'NR % 1000 == 0' | plot histogram -i - -x 9
```

<img src="figures/hist.png" width="400">

Winsorize long tails:

```
samtools view aln.bam | plot histogram -i - -x 9 -xwin 3
```

<img src="figures/hist-win.png" width="400">


### Installation and requirements

* [R](https://cran.r-project.org/) and more specifically `Rscript` which should be already installed on *.nix systems.
* R package [ggplot2](http://ggplot2.org/). To install open up the R terminal and issue: `install.packages('ggplot2')`

Once these requirements are satisfied there is no installation needed, just execute `./plot` or put the `plot` script in directory on your PATH,  *e.g.* in `~/bin/`. It might be necessary to change permission of `plot` to be executable, for this execute `chmod 755 plot`.

### Input and output

### Handling skewed data

### Handling large dataset

## TODO

* Options to label axes `-xlab <lab> -ylab <lab>`. Default to whatever appears in `-x -y`
* Option to call columns by header name
