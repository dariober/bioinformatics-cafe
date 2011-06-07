#
# Barplot of the table in Labbook 07/01/2011.
#


samfile<- c('aln21.sam', 'aln22.sam', 'aln23.sam', 'aln24.sam', 'aln25.sam', 'aln26.sam', 'aln27.sam', 'aln28.sam', 'aln29.sam', 'aln30.sam')
seqlen<- seq(21, 30)
single_map<- c(615089, 353350, 186854, 66214, 22483, 9366, 5042, 3384, 2340, 1796)
multi_map<- c(1068787, 581262, 265572, 131198, 80164, 58213, 46762, 39455, 34332, 30275)
unmapped<- c(132701, 134175, 128836, 68160, 28551, 12585, 6409, 3923, 2783, 2261)

tot_tags<- sum(single_map[1], multi_map[1], unmapped[1])

par(mar=c(5,6,4,12), cex= 0.95)
b1<- barplot(rbind(cumsum(single_map)/1000, multi_map/1000, cumsum(unmapped)/1000), 
  names.arg=seqlen, las=1, xlab= 'Tag length (bp)', col= c('dodgerblue', 'firebrick', 'darkolivegreen4'),
  yaxt='n'
  )
title(cex.main=0.95, main= 'Human TSS tags mapped to S. scrofa genome\n(macrophages monocyte derived - FANTOM5)')
axis(side= 2, at= c(0, tot_tags/1000), labels=c("",""), lwd.ticks=0)
axis(side= 2, at=seq(0 , tot_tags/1000, by=200), lwd=0, lwd.ticks=1, las= 1)
mtext(text= 'Number of tags\n(x 1000)', side= 2, line= 3)
axis(side= 4, at= c(0, tot_tags/1000), labels=c("",""), lwd.ticks=0)
axis(side= 4, at= seq(0 , tot_tags/1000, length.out=11), lwd=0, lwd.ticks=1, las= 1, labels= paste( (seq(0 , tot_tags, length.out=11)/tot_tags)*100, '%', sep=''))
legend(x= 14.5, y= par('usr')[4], legend= rev(c('Single-mappers', 'Multi-mappers', 'Unmapped')),
  fill= rev(c('dodgerblue', 'firebrick', 'darkolivegreen4')), bty= 'n', xpd= T)
savePlot('M:/Documents/LabBook/LabBook_Figures/20110110_barplot_mapping.emf', 'emf')
