#
# On transcripts assembled denovo by cufflinks
#
#

## This table created by 20100512_cuffcompare.sql
## For each dataset, has the number of exons in each transcript (no_exons) and the transcript FPKM
##
exons<- read.table('C:/Tritume/exons.txt', header=T, sep='\t')
exons[1:10,]

## Summary stats for each dataset 

exon_stats<- data.frame(dataset= unique(exons$source))
se<- function(x) sd(x)/sqrt(length(x))
funs<- c('length', 'median', 'mean', 'se')
for(fun in funs){
    exon_stats<- cbind(exon_stats,
        round(aggregate(exons$fpkm, by=list(exons$source), fun)$x, 2)
        )
    }
for(fun in funs){
    exon_stats<- cbind(exon_stats,
        round(aggregate(exons$fpkm[exons$no_exons == 1], by=list(exons$source[exons$no_exons == 1]), fun)$x, 2)
        )
    }
colnames(exon_stats)<- c('dataset', funs, paste('single', funs, sep= '_'))

exon_u<- aggregate(exons$fpkm[exons$class_code == 'u'], by=list(exons$source[exons$class_code == 'u']), length)
names(exon_u)<- c('dataset', 'n_intergenic')
exon_u1<- aggregate(exons$fpkm[exons$class_code == 'u' & exons$no_exons == 1], by=list(exons$source[exons$class_code == 'u' & exons$no_exons == 1]), length)
cbind(exon_u, n_inter_single= exon_u1$x, perc= round(exon_u1$x/exon_u$n_intergenic, 2))

aggregate(exons$fpkm[exons$class_code == 'u'], by=list(exons$source[exons$class_code == 'u']), mean)
aggregate(exons$fpkm[exons$class_code == 'u' & exons$no_exons == 1], by=list(exons$source[exons$class_code == 'u' & exons$no_exons == 1]), mean)

barplot(as.matrix(exon_stats[,c("median", "single_median")]), beside=T)

hist(log(exons$fpkm[exons$source == '20100317_RNAseq_CTRL']))




