library(RODBC)
conn_db<- odbcConnect(dsn= 'pgCAGE', case='nochange', uid= 'dberaldi', pwd= 'mypassword')

cufflinks<- sqlQuery(conn_db, "
    SELECT rname, source, feature, \"start\", \"end\", transcript_id, fpkm FROM cufflinks_transcript_gtf WHERE feature LIKE 'transcript'
    ")
odbcClose(conn_db)

## ----------------------------[ Transcriptome complexity ]-------------------------

## Do few genes contribute to most of the expression (i.e. low complexity)?
## Or no major genes make up most of transcriptome (high complexity)?

## Transcripts w/ reference
ref<- cufflinks[cufflinks$source == '20100317_RNAseq_CTRL' | cufflinks$source == '20100317_RNAseq_LPS', ]

## Plot cumulative
par(cex=0.85)
trans<- ref$fpkm[ref$source == '20100317_RNAseq_CTRL']
cumtrans<- cumsum(sort(trans, decreasing= TRUE))
plot(1:length(cumtrans), (cumtrans/sum(trans))*100, type= 'l', log='x', xlab= 'Number of genes', ylab= 'Fraction of mRNA pool (%)', 
    main= 'Transcriptome complexity', col= 'blue', cex.main= 1.25)

trans<- ref$fpkm[ref$source == '20100317_RNAseq_LPS']
cumtrans<- cumsum(sort(trans, decreasing= TRUE))
points(1:length(cumtrans), (cumtrans/sum(trans))*100, type= 'l', col= 'red')

legend('topleft', legend= c('CTRL', 'LPS'), lty= c('solid', 'solid'), col= c('blue', 'red'))
savePlot('M:/Documents/LabBook/LabBook_Figures/20100513_transcriptome_complexity.emf', 'emf')
head(cumtrans)

max(ref$fpkm[ref$source == '20100317_RNAseq_CTRL'])

# ----------------------------[ Resource allocation ]---------------------------

## This table is the one on LabBook 13 May 2010 with little re-arrangement of the
#  column names.
go<- read.table('C:/Tritume/go_terms.txt', sep='\t', header=T)
par(mar= c(4,10,2,1))
barplot(rev(go$perc_tot_expression), horiz= T, col= rep(c('lightblue', 'lightgrey', 'pink'), each= 5), 
    names.arg= substring(rev(go$GO_term), 1, 20), las=1, xlab= '', main= 'Resource allocation')
mtext(side=1, line= 2.5, text= 'Fraction of transcriptome', cex= 0.95)
text(x= 0.08, y= c(19,12,6)-3, labels= c('Biol. proc.', 'Mol. func.', 'Cell. comp.'), cex= 1.1, col=c('red', 'darkgrey', 'blue'), font=2)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100513_resource_allocation.emf', 'emf')

## -----------------------[ Tritume ]--------------------

ref<- cufflinks[cufflinks$source == '20100317_RNAseq_CTRL', ] ##  | cufflinks$source == '20100317_RNAseq_LPS'
plot(log(ref$fpkm), log((ref$end - ref$start)))
cor.test(ref$fpkm, (ref$end - ref$start))


denovo<- cufflinks[(cufflinks$source == '20100317_RNAseq_CTRL_noGTF' | cufflinks$source == '20100317_RNAseq_LPS_noGTF'), ]
denovo<- cufflinks[(cufflinks$source == '20100317_RNAseq_LPS_noGTF'), ]



dim(denovo)

ntrans<- vector(length=0, mode= 'numeric')
    rpkm<- seq(0, 100)
    for(i in rpkm){
        ntrans<- append(ntrans, length(na.omit(denovo$fpkm[denovo$fpkm > i & denovo$source == '20100317_RNAseq_CTRL_noGTF'])))
        }
    plot(rpkm, ntrans, ylim= c(0, max(ntrans)), type= 'l', main= 'Number of transcripts/RPKM (CTRL)', xlab= 'RPKM', ylab= 'N. transcripts')
    points(x= c(0,2,5,10), y= ntrans[which(rpkm %in% c(0,2,5,10))], pch= 16, col= 'red')
    text(x= c(0,2,5,10) + 5, y= ntrans[which(rpkm %in% c(0,2,5,10))], 
        labels= paste('rpkm > ', c(0,2,5,10), '\nn= ', ntrans[which(rpkm %in% c(0,2,5,10))], sep=''), 
         adj= 0, cex= 1)

