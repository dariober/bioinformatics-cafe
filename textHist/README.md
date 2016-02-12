## Text histogram

Minimalist command line tool to create histograms in ascii characters from data in files. Useful for
quickly assessing the distribution of data.

### Example

```
textHist.R -i test_data/data.txt -c 1
||||||||                                           -1.75
||||||||||                                         -1.25
||||||||||||||||||||||||||||||||||||||||           -0.75
|||||||||||||||||||||||||||||||||||||||||||||||||| -0.25
||||||||||||||||||||||||||||                       0.25
||||||||||||||||||||||||||||||||                   0.75
||||||||||||||||||                                 1.25
||||||                                             1.75
||||||||                                           2.25
+------------------------------------------------+
0                                                25
|: 2x
```

See also `run_tests.sh`

### Requires

* [R](https://cran.r-project.org/), `Rscript` in fact.
* R package [argparse](https://cran.r-project.org/web/packages/argparse/index.html)