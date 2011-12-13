#
# See labbook 26/05/2011
# Some plots relative to histone methylation
# 

library(sqldf)
library("geneplotter")  ## from BioConductor
require("RColorBrewer") ## from CRAN

# ------------------------[ Input/Output ]-------------------------------------

setwd("F:/data/20110602_markb_chipseq")
deseq<- read.table('deseq_nbinomtest.txt', header= TRUE, sep= '\t')     ## This is from deseq
annot<- read.table('con_k4b_annotation-1.txt', header= TRUE, sep= '\t', stringsAsFactors= FALSE, quote= '', comment.char= '') ## This is from Homer's annotatePeaks.pl

## Be sure deseq and annotatePeaks are using the same peak Ids
deseq_annot<- sqldf("select * from deseq inner join annot on id = PeakID")

# -----------------------------------------------------------------------------
tiff('deseq_basemean_diff.tiff',
      res= 250,                             ## Resolution in dpi
      pointsize= 6,                         ## This determines the size of the writings
      units= 'cm', width = 6, height = 5.5 ## Unit of measure and size
      )
par(las= 1)
plot( deseq_annot$baseMean, deseq_annot$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( deseq_annot$padj < .1, "red", "black" ),
     ylab= '', xlab= '', cex.main= 0.85, main= '')
abline(h= c(-2, 2), col= 'dodgerblue', lwd= 1)
mtext(side= 1, text= 'Base mean', line= 2.5)
mtext(side= 2, text= 'Log2 fold change', line= 2.5, las= 0)
mtext(side= 3, text= 'Differentially bound peaks (padj <.1)', line= 0.5)
graphics.off()


tiff('distance_to_tss_all.tiff', res= 250, pointsize= 9, units= 'cm', width = 10, height = 10)
par(las= 1)
smoothScatter(deseq_annot$Distance_to_TSS, deseq_annot$baseMean, bandwidth=0.001, main= '', pch= 19, cex= 0.5, col= 'grey60',
              xlab='', ylab='', xaxt= 'n', yaxt= 'n')
grid()
x<- seq(par('xaxp')[1], par('xaxp')[2], length.out= 9)
axis(side= 1, at= x, labels= formatC(x/1000000, format= 'f', digits= 1, big.mark= ','))
y<- seq(par('yaxp')[1], par('yaxp')[2], length.out= 6)
axis(side= 2, at= y, labels= formatC(y/1000, format= 'd', big.mark= ','))
mtext(side=1, text= 'Distance from nearest TSS (Mb)', line= 2.5)
mtext(side=2, text= 'Base mean (x1000)', line= 2.5, las= 0)
graphics.off()
## tss_dist<- sapply(deseq_annot$Distance_to_TSS, log_dist)
genes<- sqldf("select Gene_Name, Distance_to_TSS, baseMean,
               case when upper(Gene_Name) = 'CSF1R' then 'deeppink4'
                    when upper(Gene_Name) = 'NFKB2' then 'orangered4'
                    when upper(Gene_Name) = 'CSF3' then 'green4' END AS colour
               from deseq_annot where upper(Gene_Name) in ('CSF1R', 'NFKB2', 'CSF3') ")
## smoothScatter(tss_dist, deseq_annot$baseMean, bandwidth=0.001, main= '', pch= 19, cex= 0.5, col= 'grey60', xlab='Log2 TSS distance', ylab='')

tiff('distance_to_tss.tiff', res= 250, pointsize= 9, units= 'cm', width = 10, height = 10)
par(las= 1)
smoothScatter(deseq_annot$Distance_to_TSS[deseq_annot$Distance_to_TSS > -20000 & deseq_annot$Distance_to_TSS < 20000],
              deseq_annot$baseMean[deseq_annot$Distance_to_TSS > -20000 & deseq_annot$Distance_to_TSS < 20000],
              bandwidth=0.001, main= '', pch= 19, cex= 0.5, col= 'grey60',
              xlab='', ylab='', xlim= c(-20000, 20000), xaxt= 'n', yaxt= 'n')
grid()
points(genes$Distance_to_TSS, genes$baseMean, pch= 19, col= as.character(genes$colour))
text(genes$Distance_to_TSS, genes$baseMean+500, labels= genes$Gene_Name, font= 1, col= as.character(genes$colour), cex= 0.85)
x<- seq(par('xaxp')[1], par('xaxp')[2], length.out= 8)
axis(side= 1, at= x, labels= formatC(x/1000, format= 'd', big.mark= ','))
y<- seq(par('yaxp')[1], par('yaxp')[2], length.out= 6)
axis(side= 2, at= y, labels= formatC(y/1000, format= 'd', big.mark= ','))
mtext(side=1, text= 'Distance from nearest TSS (kb)', line= 2.5)
mtext(side=2, text= 'Base mean (x1000)', line= 2.5, las= 0)
## savePlot('distance_to_tss.emf', 'emf')
graphics.off()


h1<- hist(tss_dist, breaks= 100, plot= FALSE)
points(x= h1$mids, y= h1$counts, type= 'l', lwd= 2, col= 'firebrick4')

log_dist<- function(x){
      #
      if (is.na(x)){
            return(NA)
      } else if(x == 0){
            return(0)
      } else if (x < 0){
            return(-log2(abs(x)))
      } else if (x > 0){
            return(log2(x))
      } else{
            print(x)
            stop('Unexpected input')
      }
}


points(deseq_annot$Distance_to_TSS[match(c('CSF1R, NFKB2', 'CSF3'), deseq_annot$Gene_Name) ])

plot(deseq_annot$Distance_to_TSS, deseq_annot$baseMean, pch= 19)

sqldf("select * from deseq_annot where baseMean = (select max(baseMean) from deseq_annot)")
log2(deseq_annot$Distance_to_TSS[deseq_annot$Distance_to_TSS > 0])