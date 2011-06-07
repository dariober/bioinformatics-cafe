#
#  Produce a barplot of sequences aligned (total), single mappers, multi-mappers
#

## These data comes from labbook 16/12/2010
aln<- c('aln21.sam', 'aln22.sam', 'aln23.sam', 'aln24.sam', 'aln25.sam', 'aln26.sam', 'aln27.sam')
tot_reads<- c(3899341, 1858994, 1400240, 1125787, 991171, 907410, 850898)
single_mappers<- c(1926389, 380444, 190833, 90015, 58035, 44948, 41704)
multi_mappers<- c(1858994, 1400240, 1125787, 991171, 907410, 850898, 787031)
unmapped<- c(113958, 78310, 83620, 44601, 25726, 11564, 22163)

cum_single<- cumsum(single_mappers)
cum_unmapped<- cumsum(unmapped)

barplot_df<- cbind(cum_single, multi_mappers, cum_unmapped)

windows(width= 16/2.54, height= 12/2.54)
par(mar=c(5,6,4,12), cex= 0.9)

b<- barplot(t(barplot_df)/1000, names.arg='', las=1, xlab= '', col= c('dodgerblue', 'firebrick', 'darkolivegreen4'),
  yaxt='n', cex.names= 0.9, cex.axis= 0.9 )
title(cex.main=1.1, main= 'Mapping pig CAGE tags\nagainst the pig genome (S. scrofa 9)')
mtext(side= 1, text= 21:27, at= b, line= 0.5, cex= 0.9)
mtext(side= 1, text= 'Read length (bp)', line= 2, cex= 0.9)
axis(side= 2, at= c(0, tot_reads[1]/1000), labels=c("",""), lwd.ticks=0)
axis(side= 2, at=seq(0 , tot_reads[1]/1000, by=500), lwd=0, lwd.ticks=1, las= 1)
mtext(text= 'N. reads (x 1000)', side= 2, line= 3.5)
axis(side= 4, at= c(0, tot_reads[1]/1000), labels=c("",""), lwd.ticks=0)
axis(side= 4, at= seq(0 , tot_reads[1]/1000, length.out=11), lwd=0, lwd.ticks=1, las= 1, labels= paste( (seq(0 , tot_reads[1], length.out=11)/tot_reads[1])*100, ' %', sep=''))
legend('topright', inset= c(-0.7, 0),  legend= rev(c('Single-mappers', 'Multi-mappers', 'Unmapped')),
  fill= rev(c('dodgerblue', 'firebrick', 'darkolivegreen4')), bty= 'n', xpd= T)

savePlot('M:/Documents/LabBook/LabBook_Figures/20110210_cage050810_barplot_mapping.emf', 'emf')

par('usr')


