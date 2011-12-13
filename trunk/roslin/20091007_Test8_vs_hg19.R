library("BSgenome")
available.genomes()

## Download your favourite genome as a library
# NB: This takes a while (but it has to be done only once)
# source("http://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")

## Load reference genome
library(BSgenome.Hsapiens.UCSC.hg19)

## Sequences to pull out 
hg19c<- read.table("~/Documents/Tritume/dario_Test8_hg19.txt", sep="\t", 
    header=T, comment.char="")
hg19c<- hg19c[hg19c$NoMatches == 1, ]

## Extract sequences
matchSeq<- getSeq(Hsapiens, names= hg19c$chromosome,   start= hg19c$position - 1,
   end= hg19c$position + 27, strand= hg19c$strand)

   
 
## Extract the base preceeding the TSS
## This might be slow!!

tss_minus1<- function(seq, strand){
    if(strand == '+'){
        tss<- substring(seq, 1, 1)
        return(tss)
        }
    if(strand == '-'){
        tss<- substring(seq, nchar(seq), nchar(seq))
        return(tss)
        }
    }

tss_x1<- sapply(matchSeq, substring, 1, 1)

tss_x<- vector(length= length(matchSeq))
for(i in 1:length(matchSeq)){
    tss_x[i]<- tss_minus1(matchSeq[i], hg19c$strand[i])
    }   

aggregate(tss_x[hg19c$strand == '+'], list(tss_x[hg19c$strand == '+']), length)
aggregate(tss_x1[hg19c$strand == '+'], list(tss_x1[hg19c$strand == '+']), length)
