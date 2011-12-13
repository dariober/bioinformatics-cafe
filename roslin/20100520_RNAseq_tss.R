## This table created by 20100520_RNAseq_tss.sql
denovo_tss<- read.table('C:/Tritume/denovo_tss.txt', sep='\t', header=T)
par(mfrow= c(1,2))
    hist(denovo_tss$denovo_rel_position, 
        main= '',
        xlab= 'de novo relative position\n(negative: de novo starts earlier)')
    legend('topleft', legend= paste('n=', length(denovo_tss$denovo_rel_position)), bty= 'n')
    hist(denovo_tss$denovo_rel_position[denovo_tss$denovo_rel_position < 100 & denovo_tss$denovo_rel_position > -100], 
        breaks= 20, main=NA, xlab= 'Position between -100 and +100\n(negative: de novo starts earlier)')
    mtext(side= 3, text='Position of 5\' ends of de novo transcripts relative to reference transcripts', outer=T, line= -2, cex= 1.25)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100520_5prime_end.emf', 'emf')