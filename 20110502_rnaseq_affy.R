library(RODBC)
conn<- odbcConnect(dsn= 'pgVitelleschi')

platforms<- sqlQuery(conn, 'select * from platforms')
platforms[1:10,]

bmp('M:/Documents/LabBook/LabBook_Figures/20110502_rnaseq_affy.bmp', height= 7, width= 17, units= 'cm', res= 150, pointsize = 12)
par(pch= 19, cex.axis= 0.85, col= 'grey25', las= 1, mfrow= c(1,3), mar=c(4,2.5,3,0), oma= c(2,2,0,0.1), mgp=c(2,1,0))
    for(i in c(20,5,2)){
    plot(platforms$affy_log_fc, platforms$rnaseq_log_fc,
        xlab= 'Affymetrix log2', ylab= '',
        type= 'n',
        xlim= c(-i,i), ylim= c(-i,i)
    )
    grid()
    abline(a= 0, b= 1, col= 'steelblue', lwd=2)
    points(platforms$affy_log_fc, platforms$rnaseq_log_fc, cex= 0.2)
    lines(loess.smooth(platforms$affy_log_fc, platforms$rnaseq_log_fc), col= 'red', lwd= 2)
}
mtext(text= 'RNAseq log2', outer= TRUE, side= 2, las= 0, cex= 0.66)
dev.off()


## Significant in Affy and/or RNAseq (cuffdiff)
colv<- ifelse(platforms$avg_fdr < 0.01 & platforms$significant == 'yes', 'red',                          ## Both sign.
           ifelse(platforms$avg_fdr < 0.01 & platforms$significant == 'no', 'blue',                      ## affy only
               ifelse(platforms$avg_fdr > 0.01 & platforms$significant == 'yes', 'green', 'transparent')      ## RNAseq only
           ) 
       )
ncounts<- aggregate(colv, by= list(colv), length)
leg<- paste(c('Both   ',
              'Affy   ',
              'RNAseq ',
              'NS     '),
            'n= ', c(ncounts$x[ncounts$Group.1 == 'red'],
                     ncounts$x[ncounts$Group.1 == 'blue'],
                     ncounts$x[ncounts$Group.1 == 'green'],
                     ncounts$x[ncounts$Group.1 == 'transparent']),
            sep= '')
bmp('M:/Documents/LabBook/LabBook_Figures/20110502_rnaseq_affy_sign.bmp', height= 9, width= 9, units= 'cm', res= 300, pointsize = 8)
par(pch= 21, cex.axis= 0.85, col= 'grey60', las= 1, mgp=c(2,1,0))
plot(platforms$affy_log_fc, platforms$rnaseq_log_fc,
        xlab= 'Affymetrix log2', ylab= 'RNAseq log2', main= 'Transcripts found D.E. by Affymetrix and RNAseq', cex.main= 0.85,
        type= 'n',
        xlim= c(-8,8), ylim= c(-8,8)
    )
grid()
points(platforms$affy_log_fc, platforms$rnaseq_log_fc, cex= 0.5, col= colv, pch= 19)
par(font= 21)
legend('topleft', legend= leg, box.col= 'transparent',
    pch= 19, col= c('red', 'blue', 'green', 'transparent'), text.col= 'black',
    cex= 0.8, bg= 'white')
box()
dev.off()

plot(platforms$rnaseq_ctrl, platforms$rnaseq_lps, col= colv, pch= 19, log= 'xy', cex= 0.5)
grid()
par(font= 24)
legend('topleft', legend= leg, box.col= 'transparent',
    pch= 19, col= c('red', 'blue', 'green', 'transparent'), text.col= 'black',
    cex= 0.8, bg= 'white')
