#
# CD14 expression - RNASeq + CAGE
#


library(rpileup)
library(RODBC)
conn<- odbcConnect(dsn= 'pgVitelleschi')
cage<- sqlQuery(conn, "
    select * from cage050810_tss where rname = '2' and tss_pos between 129356198-500 and 129357848+500;
")

bams<- c(
    'F:/data/bamfiles/20110202_rnaseq_am_ctrl/accepted_hits.bam',
    'F:/data/bamfiles/20110202_rnaseq_am_lps/accepted_hits.bam',
    'F:/data/bamfiles/20110202_rnaseq_bmdm_ctrl/accepted_hits.bam',
    'F:/data/bamfiles/20110202_rnaseq_bmdm_lps/accepted_hits.bam'
)
piles<- import.pileup(bams, rname= 2, from= 129356198-500, to= 129357848+500)
## plot.pileup(piles[piles$pos > 129357500 & piles$pos < 129358000, ], overplot= F)

dev.off()
windows(width= 16/2.54, height= 12/2.54)
plot.pileup(piles, overplot= F, col= c('steelblue', 'salmon'), type= 'h', cex= 0.8,
            pileup.names= c('am - ctrl', 'am - lps', 'bmdm - ctrl', 'bmdm - lps'), oma= c(4,4,4,1))

############
library(biomaRt)
ensembl<- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
mart<- useDataset("sscrofa_gene_ensembl", ensembl)
writeClipboard(listAttributes(mart)[,1])
writeClipboard(listFilters(mart)[,1])

## Get transcript features in a genomic window:
features<- getBM(mart= mart,
    attributes= c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_exon_id', 'transcript_start', 'transcript_end', 'strand' ,'exon_chrom_start', 'exon_chrom_end', '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end'),
    filters= c('chromosome_name', 'start', 'end'),
    value= list(2, par('usr')[1], par('usr')[2])
)

#exons<- c(129356198, 129357428, 129357509, 129357848)
#for(i in 1:4){
#    par(mfg=c(i, 1))
#    abline(v= exons, 
#       col= 'blue', lty= 'dotted')
#}

par(mfg=c(1,1))
nlines<- 3 ## distance from plot border to start drawing transcripts

## Plot entire transcript as a line
abline(v=1, xpd= TRUE, col= 'transparent'); abline(v=1, xpd= NA, col= 'transparent')
segments(x0= min(features$transcript_start), x1= max(features$transcript_end),
         y0= par('usr')[4] + (strheight('x')*nlines), y1= par('usr')[4] + (strheight('x')*nlines),
         lwd= 2, col= 'firebrick4', xpd= NA
)


box.cex<- 1.05
## Plot exon blocks
rect(xleft= features$exon_chrom_start,
     xright= features$exon_chrom_end,
     ybottom= (par('usr')[4] + (strheight('x')*nlines)) * (1/box.cex),
     ytop= (par('usr')[4] + (strheight('x')*nlines)) * box.cex,
     border= 'firebrick4', col= 'firebrick4', xpd= NA
)

## Plot 5'-UTR blocks
rect(xleft= features$"5_utr_start",
     xright= features$"5_utr_end",
     ybottom= (par('usr')[4] + (strheight('x')*nlines)) * (1/box.cex),
     ytop= (par('usr')[4] + (strheight('x')*nlines)) * box.cex,
     border= 'firebrick4', col= 'white', xpd= NA
)

## Plot 3'-UTR blocks   
rect(xleft= features$"3_utr_start",
     xright= features$"3_utr_end",
     ybottom= (par('usr')[4] + (strheight('x')*nlines)) * (1/box.cex),
     ytop= (par('usr')[4] + (strheight('x')*nlines)) * box.cex,
     border= 'firebrick4', col= 'white', xpd= NA
)

## CAGE peaks
for(i in 1:4){
    par(mfg=c(i, 1))
    points(x= cage$tss_pos, y= cage$tss_count*1.5, type= 'h', col= 'black', lwd=1.8)
}

## Affymetrix probes
affy_probes<- getBM(mart= mart,
    attributes= 'affy_porcine',
    filters= c('chromosome_name', 'start', 'end'),
    value= list(2, par('usr')[1], par('usr')[2])
)

title(main= 'CD14 - RNAseq and CAGE expression', outer= TRUE)
savePlot('M:/Documents/LabBook/LabBook_Figures/20110413_CD14_expression.emf', 'emf')
savePlot('M:/Documents/LabBook/LabBook_Figures/20110413_CD14_expression.bmp', 'bmp')
graphics.off()
names(cage)