#
# See labbook 26/05/2011
# Some plots relative to histone methylation
# 

library(sqldf)
library("geneplotter")  ## from BioConductor
require("RColorBrewer") ## from CRAN
library(biomaRt)

# ------------------------[ Input/Output ]-------------------------------------

setwd("F:/data/20110602_markb_chipseq")
deseq<- read.table('deseq_nbinomtest.txt', header= TRUE, sep= '\t')     ## This is from deseq
annot<- read.table('con_k4b_annotation-1.txt', header= TRUE, sep= '\t', stringsAsFactors= FALSE, quote= '', comment.char= '') ## This is from Homer's annotatePeaks.pl

# -----------------------------------------------------------------------------
## Be sure deseq and annotatePeaks are using the same peak Ids
deseq_annot<- sqldf("select * from deseq inner join annot on id = PeakID")

sqldf("select count(distinct Gene_Name)from deseq_annot where padj < 0.1 and Gene_Name not like '' ") ## Number of genes (gene names) D.E.= 165
sqldf("select count(distinct PeakID)from deseq_annot where padj < 0.1 ") ## Number of Peaks D.E.


## Differential bound peaks
tiff('deseq_basemean_diff.tiff', res= 250, pointsize= 6, units= 'cm', width = 6, height = 5.5)
par(las= 1)
plot( deseq_annot$baseMean, deseq_annot$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( deseq_annot$padj < .1, "red", "black" ),
     ylab= '', xlab= '', cex.main= 0.85, main= '')
abline(h= c(-2, 2), col= 'dodgerblue', lwd= 1)
mtext(side= 1, text= 'Base mean', line= 2.5)
mtext(side= 2, text= 'Log2 fold change', line= 2.5, las= 0)
mtext(side= 3, text= 'Differentially bound peaks (padj <.1)', line= 0.5)
graphics.off()

## All peaks distance to TSS 
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
genes<- sqldf("select Gene_Name, Distance_to_TSS, baseMean,
               case when upper(Gene_Name) = 'CSF1R' then 'deeppink4'
                    when upper(Gene_Name) = 'NFKB2' then 'orangered4'
                    when upper(Gene_Name) = 'CSF3' then 'green4' END AS colour
               from deseq_annot where upper(Gene_Name) in ('CSF1R', 'NFKB2', 'CSF3') ")

## Peaks within -20 and +20 kb from nearest TSS.
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

# ----------------------------[ BiomaRt GO terms ]-----------------------------

ensembl<- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
mart<- useDataset("mmusculus_gene_ensembl", ensembl)
sig_peaks<- sqldf(" select * from deseq_annot where padj < 0.1 ")
go_terms<- getBM(
             attributes= c("entrezgene", "external_gene_id", "name_1006", "goslim_goa_accession", "goslim_goa_description"), 
             filters= "entrezgene",
             value= sig_peaks$Entrez_ID,
             mart= mart)
genes_to_go<- sqldf("select distinct entrezgene, external_gene_id as gene_name, name_1006 as go_term_bp from go_terms order by go_term_bp, gene_name; ")
write.table(genes_to_go, 'genes_to_go.txt', row.names= FALSE, col.names= TRUE, sep= '\t')

sqldf(" select count(distinct goslim_goa_accession), count(distinct name_1006), count(distinct entrezgene) from go_terms ")
no_genes<- sqldf(" select count(distinct Entrez_ID) from sig_peaks ")[1,1]
go_count<- sqldf(" select name_1006 as go_term_bp, count(distinct entrezgene) as no_distinct_genes from go_terms group by name_1006 order by no_distinct_genes desc")
write.table(go_count, 'go_terms.txt',
            row.names= FALSE, col.names= TRUE, sep= '\t')
graphics.off()
windows(height= 9/2.54, width= 14/2.54)
par(mar= c(4,9,4,1), cex= 0.85, mgp= c(3, 0.6, 0))
barplot(go_count$no_distinct_genes[1:10], horiz= TRUE, names.arg= go_count$go_term_bp[1:10], las= 1,
        main= '', col= 'salmon')
axis(side= 3, at= seq(0, 0.5, 0.1)*no_genes, labels=  paste((seq(0, 0.5, 0.1)*100), '%', sep= ''))
mtext(side= 1, line= 2, text= 'No. genes', cex= 0.85)
mtext(side= 3, line= 2.75, text= 'Genes near differential peaks: Associated GO terms', font= 2, cex= 0.85)
savePlot('go_terms.emf', 'emf')

# -------------------------[ TRITUME ]-----------------------------------------
goslim_count<- sqldf(" select goslim_goa_description as goslim_descritpion, count(distinct entrezgene) as no_distinct_genes from go_terms group by goslim_goa_description order by no_distinct_genes desc")

write.table(goslim_count, 'goslim_terms.txt',
            row.names= FALSE, col.names= TRUE, sep= '\t')

write.table(deseq_annot, 'deseq_annot.txt', sep= '\t', col.names= TRUE, row.names= FALSE)

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


