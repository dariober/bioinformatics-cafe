Fast extraction of sequences from FASTA file
============================================

`fastaFromBed.py` accomplishes the same task as bedtools [fastaFromBed](http://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html)
but it is much faster as the reference fasta sequences is loaded in memory. This of course comes at the expense 
of spending time loading the reference and keeping it in memory. 
The speed up is most noticeable when millions of sequences needs to be extracted. Approximately,
It takes about 1 min and 4GB of memory to extract 1 million sequences from the mouse genome mm9 (2.6GB).

Help (*might be outdated*):

```
 ./fastaFromBed.py -h
usage: fastaFromBed.py [-h] [--fasta FASTA] [--bed BED] --fo FO [--name]
                       [--tab] [--strand] [--verbose] [--version]

DESCRIPTION
    Extract DNA sequences into a fasta file based on feature coordinates.
    
    Same jobs as bedtools fastaFromBed with the difference that the fasta file is loaded
    in memory. This makes fastaFromBed.py much faster than bedtools when extracting
    several sequences (e.g. millions).
    
    Note: It takes about 1 min and 4GB of memory to extract 1 million sequences from
    the mouse genome mm9 (2.6GB)

    See also 
   

optional arguments:
  -h, --help            show this help message and exit
  --fasta FASTA, -fi FASTA
                        Input fasta file. Can be gzipped (but slower)
                        
                                           
  --bed BED, -bed BED   BED file of ranges to extract from -fi. Use - to read
                        from stdin.
                        
                                           
  --fo FO, -fo FO       Output file (can be FASTA or TAB-delimited). Use -
                        for stdout.
                                            
  --name, -name         Use the name field for the FASTA header.
                                
                                           
  --tab, -tab           Write output in TAB delimited format. Default is FASTA format.
                               
                                           
  --strand, -s          Force strandedness. If the feature occupies the antisense,
                        strand, the sequence will be reverse complemented.
                        - By default, strand information is ignored.
                                           
  --verbose, -V         Show some progress report
                                           
  --version             show program's version number and exit
```