# -----------------------------------------------------------------------------
# Extract TSS sequence given coordinates in BED file(s)
# -----------------------------------------------------------------------------

library(BSgenome)

## Source genome
library("BSgenome.Hsapiens.UCSC.hg19")
cur_genome<- Hsapiens

## BED file
bed_file<- 'F:/data/20110105_FANTOM5/Macrophage%20-%20monocyte%20derived%20donor1.CNhs10861.11232-116C8.hg19.ctss.bed'


## Import data
bed<- read.table(file= bed_file, sep='\t', stringsAsFactors= FALSE)
colnames(bed)<- c('rname', 'start', 'end', 'feature_id', 'score', 'strand')

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

head(bed)
120 240 000
write.table(bed, file= 'D:/Tritume/getseq.txt', sep= '\t')
object.size(bed)
is.numeric(bed$start)
object.size(tss)
length(tss)

smallbed<- bed[0, ]
for(x in unique(bed$rname) ){
    smallbed<- rbind(smallbed, head(bed[bed$rname == x, ]) )
    smallbed<- rbind(smallbed, tail(bed[bed$rname == x, ]) )
    }
dim(smallbed)

xnames<- rep(c('a', 'b', 'c'), times= c(2,4,5))
clen<- c(10, 20, 30)
names(clen)<- c('a', 'b', 'c')

merge(xnames, clen)

xlen<- clen[which(names(clen) == xnames) ]
?match
object.size(bed) ## 145,321,224 bytes
z<- import.bed(bed_file)
object.size(z) ## 
head(z)
z$space[1:10]

head(bed)

cur_genome

##WRONG:
tss<- getSeq(x= Hsapiens,
             names= tail(df$V1, n= 1000),
             start= tail(ifelse(df$V6 == '+', (df$V2 + 1), ( (df$V2+1) - 29) ), n= 1000),
             end= tail(ifelse(df$V6 == '+', (df$V2 + 1 + 29), (df$V2 + 1) ), n= 1000),
             strand= tail(df$V6, n= 1000)
             )


write.table(file= 'D:/Tritume/check_seq.txt', cbind(tail(df, n= 7),
                                                    seq= tss,
                                                    f_start= tail(ifelse(df$V6 == '+', (df$V2 + 1), (df$V2 - 29)), n= 7),
                                                    f_end= tail(ifelse(df$V6 == '+', (df$V2 + 1 + 29), (df$V2 + 1) ), n= 7),
                                                    strand= tail(df$V6, n= 7)
                                                    ),
            quote= FALSE, row.names=FALSE, sep= '\t'
            )

tss<- getSeq(x= Hsapiens,
             names= head(df$V1, n= 7),
             start= head(ifelse(df$V6 == '+', (df$V2 + 1), (df$V2 - 29)), n= 7),
             end= head(ifelse(df$V6 == '+', (df$V2 + 1 + 29), (df$V2 + 1) ), n= 7),
             strand= head(df$V6, n= 7)
             )
write.table(file= 'D:/Tritume/check_seq.txt', cbind(head(df, n= 7),
                                                    seq= tss,
                                                    f_start= head(ifelse(df$V6 == '+', (df$V2 + 1), (df$V2 - 29)), n= 7),
                                                    f_end= head(ifelse(df$V6 == '+', (df$V2 + 1 + 29), (df$V2 + 1) ), n= 7),
                                                    strand= head(df$V6, n= 7)
                                                    ),
            quote= FALSE, row.names=FALSE, sep= '\t'
            )

head(df)


getSeq(x= Hsapiens,
             names= 'chrY',
             start= 59361889 - 29, ## = 59361860
             end= 59361889 + 1,    ## = 59361890
             strand= '-'
             )

R
59362113
59362142
bed
         V1       V2       V3                        V4 V5 V6
1816572 chrY 59362141 59362142 chrY:59362141..59362142,-  1  -

