<!-- MarkdownTOC -->

- [Description](#description)
- [Examples](#examples)
- [Command line options](#command-line-options)
- [Download, Installation, Requirements](#download-installation-requirements)

<!-- /MarkdownTOC -->

Description
===========

Search a fasta file for matches to a regex and return a bed file with the
coordinates of the match and the matched sequence itself.

By default, `fastaRegexFinder.py` searches for putative G-quadruplexes on
forward and reverse strand using the quadruplex rule described at 
[Quadruplex prediction techniques](http://en.wikipedia.org/wiki/G-quadruplex#Quadruplex_prediction_techniques).
The default regex is `'([gG]{3,}\w{1,7}){3,}[gG]{3,}'`.

Output bed file has columns:

```
1. Name of fasta sequence (e.g. chromosome)
2. Start of the match
3. End of the match
4. ID of the match
5. Length of the match
6. Strand 
7. Matched sequence as it appears on the forward strand
```

For matches on the reverse strand it is reported the start and end position on
the forward strand and the matched string on the forward strand (so the G4
`GGGAGGGT` present on the reverse strand is reported as `ACCCTCCC`).

Note: Fasta sequences (chroms) are read in memory one at a time along with the
matches for that chromosome. The order of the output is: chroms as they are
found in the inut fasta, matches sorted within chroms by positions.

Examples
========

```
## Test data:
echo '>mychr' > /tmp/mychr.fa
echo 'ACTGnACTGnACTGnTGAC' >> /tmp/mychr.fa

fastaRegexFinder.py -f /tmp/mychr.fa -r 'ACTG'
    mychr   0   4   mychr_0_4_for   4   +   ACTG
    mychr   5   9   mychr_5_9_for   4   +   ACTG
    mychr   10  14  mychr_10_14_for 4   +   ACTG

fastaRegexFinder.py -f /tmp/mychr.fa -r 'ACTG' --maxstr 3
    mychr   0   4   mychr_0_4_for   4   +   ACT[3,4]
    mychr   5   9   mychr_5_9_for   4   +   ACT[3,4]
    mychr   10  14  mychr_10_14_for 4   +   ACT[3,4]

less /tmp/mychr.fa | fastaRegexFinder.py -f - -r 'A\w\wGn'
    mychr   0   5   mychr_0_5_for   5   +   ACTGn
    mychr   5   10  mychr_5_10_for  5   +   ACTGn
    mychr   10  15  mychr_10_15_for 5   +   ACTGn
```

See also this [Help simulating RRBS reads](https://www.biostars.org/p/143845/) to 
extract RRBS fragments.

Command line options
====================

*Might be outdated*

```
-h, --help            show this help message and exit
--fasta FASTA, -f FASTA
                    Input fasta file to search. Use '-' to read the file from stdin.
                                                       
                                       
--regex REGEX, -r REGEX
                    Regex to be searched in the fasta input.
                    Matches to the reverse complement will have - strand.
                    The default regex is '([gG]{3,}\w{1,7}){3,}[gG]{3,}' which searches
                    for G-quadruplexes.                                   
                                       
--matchcase, -m       Match case while searching for matches. Default is
                    to ignore case (I.e. 'ACTG' will match 'actg').
                                       
--noreverse           Do not search the reverse complement of the input fasta.
                    Use this flag to search protein sequences.                                   
                                       
--maxstr MAXSTR       Maximum length of the match to report in the 7th column of the output.
                    Default is to report up to 10000nt.
                    Truncated matches are reported as <ACTG...ACTG>[<maxstr>,<tot length>]
                                       
--seqnames SEQNAMES [SEQNAMES ...], -s SEQNAMES [SEQNAMES ...]
                    List of fasta sequences in --fasta to
                    search. E.g. use --seqnames chr1 chr2 chrM to search only these crhomosomes.
                    Default is to search all the sequences in input.
                                       
--quiet, -q           Do not print progress report (i.e. sequence names as they are scanned).                                   
                                       
--version, -v         show program's version number and exit
```

Download, Installation, Requirements
====================================

Requires only Python 2.7. To download and make executable you can use:

```
wget https://github.com/dariober/bioinformatics-cafe/blob/master/fastaRegexFinder/fastaRegexFinder.py?raw=true -O fastaRegexFinder.py
chmod a+x fastaRegexFinder.py
```

There is no installation needed, but for convenience you can copy ` fastaRegexFinder.py`
to a directory on your PATH, *e.g.* `/usr/local/bin/` or `~/bin/`.


TODO
====

* Use [https://github.com/mdshw5/pyfaidx](pyfaidx) to retrieve sequences



