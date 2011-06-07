# -----------------------------------------------------------------------------
# Extract TSS sequence given coordinates in BED file(s)
# -----------------------------------------------------------------------------

library(BSgenome)

## Source genome
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Ggallus.UCSC.galGal3")

cur_genome<- Ggallus ## Hsapiens

## BED file

## bed_file<- 'F:/data/20110105_FANTOM5/Macrophage%20-%20monocyte%20derived%20donor1.CNhs10861.11232-116C8.hg19.ctss.bed'
bed_file<- 'F:/UPDATE_010/f5pipeline/chicken.primary_cell.hCAGE/Chicken%20Aortic%20Smooth%20Muscle%20cells%20donor1.CNhs11296.11298-117B2.galGal3.ctss.bed'

## Import data
bed<- read.table(file= bed_file, sep='\t', stringsAsFactors= FALSE)
colnames(bed)<- c('rname', 'start', 'end', 'feature_id', 'score', 'strand')

## dim(bed)
## bed_bk<- bed
## --------------------------------[ Delete ]----------------------------------
#smallbed<- bed[0, ]
#for(x in unique(bed$rname) ){
#    smallbed<- rbind(smallbed, head(bed[bed$rname == x, ]) )
#    smallbed<- rbind(smallbed, tail(bed[bed$rname == x, ]) )
#    }
#bed<- smallbed
#xv<- data.frame('chrM', 16561, 16562, 'chrM_test', 1, '+')
#colnames(xv)<- colnames(bed)
#bed<- rbind(bed, xv) ## TSS sequence coordinate will be beyond chr boundary (len chrM= 16571)

## --------------------------------[ end of delete ]----------------------------

## Chromosome lengths
chrlen<- seqlengths(cur_genome)
chrlen<- data.frame(rname= names(chrlen), len= chrlen)

## Prepare genomic coordinates:
seqstart<- ifelse(bed$strand == '+', (bed$start + 2), (bed$start - 29))
seqstart[which(seqstart < 1)]<- 1

bed$seqstart<- seqstart
rm(seqstart)

seqend<- ifelse(bed$strand == '+', (bed$start + 31), (bed$start) )
bed$seqend<- seqend
rm(seqend)

bed<- merge(bed, chrlen, by= intersect("rname", "rname"), sort=FALSE)

## bed now looks like this:
#    rname     start       end                  feature_id score strand  seqstart    seqend       len
# 1   chr1 148103857 148103858 chr1:148103857..148103858,-     1      - 148103828 148103857 200994015
# 2   chr1 148109605 148109606 chr1:148109605..148109606,-     1      - 148109576 148109605 200994015
# 3   chr1 148110268 148110269 chr1:148110268..148110269,-     1      - 148110239 148110268 200994015
# 4   chr1 148103686 148103687 chr1:148103686..148103687,-     1      - 148103657 148103686 200994015

## Remove lines where the TSS is past the end of the chromosome (this is a bug in
## RIKEN's pipeline)
nrow_prefilter<- nrow(bed)
bed<- bed[which(bed$seqstart < bed$len), ]
nrow_postfilter<- nrow(bed)

## When the end point of the sequence to retrieve is > then chr. length, reset it
## to the end of the chr itself.
bed$seqend[which(bed$seqend > bed$len)]<- bed$len[which(bed$seqend > bed$len)]

# -----------------------------------------------------------------------------
# Get sequences
# -----------------------------------------------------------------------------

tss<- getSeq(x= cur_genome,
             names= bed$rname,
             start= bed$seqstart,  ## ifelse(bed$strand == '+', (bed$start + 2), (bed$start - 29)),
             end=   bed$seqend,    ## ifelse(bed$strand == '+', (bed$start + 31), (bed$start) ),
             strand=bed$strand     ## bed$strand
             )

# -----------------------------------------------------------------------------
#                                 TRITUME
# -----------------------------------------------------------------------------

getseq_file<- file('D:/Tritume/getseq.txt')
open(getseq_file)
sink("D:/Tritume/getseq.txt")
1:10
sink()
close(getseq_file)

sink("D:/Tritume/sink-examp.txt")
i <- 1:10
outer(i, i, "*")
sink()

capture.output(1:10, "D:/Tritume/sink-examp.txt")


system.time({
tss<- getSeq(x= cur_genome,
             names= bed$rname[1:500000],
             start= bed$seqstart[1:500000],  ## ifelse(bed$strand == '+', (bed$start + 2), (bed$start - 29)),
             end=   bed$seqend[1:500000],    ## ifelse(bed$strand == '+', (bed$start + 31), (bed$start) ),
             strand=bed$strand[1:500000],     ## bed$strand
             as.character= TRUE
             )
})

