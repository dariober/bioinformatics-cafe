library(edgeR)

setwd('F:/data/20100826')

## These objects created by 20100413_affyarray_RNAseq.R and 
## 20100408_RNAseq_edgeR_transcript_expr.R
load('detags500.com')
load('d')
load('de_avg')
load('lm_de2')

# save(detags500.com, file= 'detags500.com')
# save(d, file= 'd')
# save(de_avg, file= 'de_avg')
# save(lm_de2, file= 'lm_de2')

par(mfrow= c(2,1), cex= 0.8, mar= c(4,4,2,1))
    ## edgeR plot
    plotSmear(d, de.tags = detags500.com, main = "Differential expression", xlim= c(-25, -5), ylim= c(-7,7))
       abline(h = c(-2, 2), col = "dodgerblue", lwd = 2, xlab= '')
    rect(xpd= TRUE, xleft= -20, ybottom= -11.5, xright= -10, ytop= -9.5, col= 'white',
        border= 'white')
    mtext(side= 1, line= 2, text= 'logConc', cex= 0.8)
    nde<- length(detags500.com)
    n<- nrow(d$counts)
    legend('bottomleft', legend= c(
        paste('All transcripts: n= ', n), 
        paste('Transcripts with p<0.01 for d.e.: n= ', nde)
        ),
        pch= 19, col= c('black', 'red'), bty= 'n', bg= 'white', box.col= 'white',
     y.intersp= 1.5)

    ## Array vs RNAseq
    plot(de_avg$rnaseq_log2fc, de_avg$affy_log2fc, 
      ylim= c(min(de_avg$rnaseq_log2fc, de_avg$affy_log2fc), max(de_avg$rnaseq_log2fc, de_avg$affy_log2fc)),
      xlim= c(min(de_avg$rnaseq_log2fc, de_avg$affy_log2fc), max(de_avg$rnaseq_log2fc, de_avg$affy_log2fc)),
      xlab= '', ylab= 'logFC: pig Affymetrix', main= 'RNAseq vs Affymetrix', cex.main= 1.25,
      cex= 0.95)
    mtext(side= 1, line= 2, text= 'logFC: RNAseq', cex= 0.8)
    points(de_avg$rnaseq_log2fc[de_avg$affy_fdr>0.05], de_avg$affy_log2fc[de_avg$affy_fdr>0.05], col='blue')
    points(de_avg$rnaseq_log2fc, lm_de2$fitted.values, type='l', col= 'red', lwd= 2, lty= 'dotted')
    ##
    ## This points() highlights the genes from the previuos plot
    ##
    ## points(de_avg$rnaseq_log2fc[de_avg$transcript_id %in% sel_t], 
    ##     de_avg$affy_log2fc[de_avg$transcript_id %in% sel_t], col= 'red', pch= 16, cex= 0.75)
    abline(a= 0, b= 1, lty= 'dotted')
    legend('topleft', bty='n', legend= paste('Affy probes with\np > 0.05 (n= ', length(de_avg$affy_log2fc[de_avg$affy_fdr>=0.05]), ')', sep=''),
        pch= 16, col= 'blue')
    legend('bottomright', bty='n', legend= paste('n transcripts=', nrow(de_avg)))
    grid()
savePlot('rnaseq_vs_affy.emf', 'emf')
savePlot('rnaseq_vs_affy.jpeg', 'jpeg')