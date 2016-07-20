#!/bin/bash

docstring="\n
Interleave a pair of fastq files.\n
USAGE\n
interleaveFastq.sh r1.fq r2.fq\n
\n
Input can be gzipped (ending .gz) or bzip2'd (endong .bz2).
"

if [[ "$1" == "-h" ]]
then
    echo -e $docstring
    exit 1
fi

xcat1='cat '$1
if [[ "$1" == *.gz ]]
then
    xcat1='gunzip -c '$1
elif [[ "$1" == *.bz2 ]]; then
    xcat1='bunzip2 -c '$1
fi

xcat2='cat '$2
if [[ "$2" == *.gz ]]
then
    xcat2='gunzip -c '$2  
elif [[ "$2" == *.bz2 ]]; then
    xcat2='bunzip2 -c '$2
fi

paste <($xcat1 | paste - - - -) \
      <($xcat2 | paste - - - -) \
| tr '\t' '\n'