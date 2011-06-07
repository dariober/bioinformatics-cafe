#
#  Boxplots of read quality for RNAseq_20100614
#  The input ptoduced by fastx (see /exports/work/vet_roslin_nextgen/dario/fastq/20100614_RNAseq_pig_BMDM/478)
#

setwd('D:/Tritume')
ctrl_1<- read.table('BMDM_CTRL_14JUN10_300BP_s_7_1_sequence.quality.csv', header=T)
ctrl_2<- read.table('BMDM_CTRL_14JUN10_300BP_s_7_2_sequence.quality.csv', header=T)
lps_1<- read.table('BMDM_LPS_14JUN10_300BP_s_8_1_sequence.quality.csv', header=T)
lps_2<- read.table('BMDM_LPS_14JUN10_300BP_s_8_2_sequence.quality.csv', header=T)

par(mfrow= c(2,2), mar=c(1,1,0,0), oma= c(4,4,4,3), cex= 0.85)
  bxp(list(stats= t(ctrl_1[,c('lW', 'Q1', 'med', 'Q3', 'rW')]), n= ctrl_1$count), xaxt='n', ylim=c(0,40), frame.plot=FALSE, boxfill=c(rep('transparent',4), 'lightblue'))
  axis(side=3, at=seq(0, 57, by= 5), labels= c(1, seq(5, 57, by= 5) ) )
  abline(h=c(10, 15, 20), lty= 'dotted', col= c('blue', 'red', 'blue'))
  legend('bottomleft', legend='CTRL-1', bty='n', text.col= 'darkblue', cex=1.2)

  bxp(list(stats= t(ctrl_2[,c('lW', 'Q1', 'med', 'Q3', 'rW')]), n= ctrl_2$count), yaxt='n', xaxt='n', ylim=c(0,40), frame.plot=FALSE, boxfill=c(rep('transparent',4), 'lightblue'))
  axis(side=3, at=seq(0, 57, by= 5), labels= c(1, seq(5, 57, by= 5) ) )
  axis(side=4, at=seq(0, 40, by= 10))
  abline(h=c(10, 15, 20), lty= 'dotted', col= c('blue', 'red', 'blue'))
  legend('bottomleft', legend='CTRL-2', bty='n', text.col= 'darkblue', cex=1.2)
  
  bxp(list(stats= t(lps_1[,c('lW', 'Q1', 'med', 'Q3', 'rW')]), n= lps_1$count), ylim=c(0,40), frame.plot=FALSE, xaxt='n', boxfill=c(rep('transparent',4), 'lightblue'))
  axis(side=1, at=seq(0, 57, by= 5), labels= c(1, seq(5, 57, by= 5) ) )
  abline(h=c(10, 15, 20), lty= 'dotted', col= c('blue', 'red', 'blue'))
  mtext(side=1, text= 'read position', line= 2.5)
  mtext(side=2, text= 'Quality score (Solexa)',  line= 2.5)
  legend('bottomleft', legend='LPS-1', bty='n', text.col= 'firebrick4', cex=1.2)

  bxp(list(stats= t(lps_2[,c('lW', 'Q1', 'med', 'Q3', 'rW')]), n= lps_2$count), yaxt='n', ylim=c(0,40), frame.plot=FALSE, xaxt='n', boxfill=c(rep('transparent',4), 'lightblue'))
  axis(side=1, at=seq(0, 57, by= 5), labels= c(1, seq(5, 57, by= 5) ) )
  axis(side=4, at=seq(0, 40, by= 10))
  abline(h=c(10, 15, 20), lty= 'dotted', col= c('blue', 'red', 'blue'))
  legend('bottomleft', legend='LPS-2', bty='n', text.col= 'firebrick4', cex=1.2)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100113_rnaseq_20100614_boxplot.emf', 'emf')

