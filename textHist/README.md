## Text histogram

Command line tool to create histograms from data files. Output is in ascii characters. Useful for
quickly scanning data.

### Example

```
 ./textHist.R -i test_data/data.txt -c 1
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
```

See also `run_tests.sh`

### Requires

* R
* R package `argparse`