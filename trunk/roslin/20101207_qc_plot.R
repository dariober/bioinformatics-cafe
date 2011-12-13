#
#  Quality charts for CAGE library 'CAGE 05/08/2010' sequence returned on 6/12/10
# --------------------------------------------------------------------------------

library(RODBC)

# sumstats is the output of 
# fastx_quality_stats \
#    -i /exports/work/vet_roslin_nextgen/dario/fastq/20100805_CAGE_AM/Beraldi_CAGE_50810_101202_EBRI093151_0060_s_8_sequence.txt \
#    -o 20101206_CAGE_qc.txt

sumstats<- read.table('F:/data/20101206_CAGE/20101206_CAGE_qc.txt', sep= '\t', header= T)
head(sumstats)

phred2p<- function(Q, d){
    # Returns the p value associated with a phred score
    p<- formatC(10^(-Q/10), digits= 10, format= 'f')
    p<- substr(p, 2, d+2)
    return(p)
    }

par(mfrow= c(3,1) , mar=c(0,5,1,5), oma= c(4,0,2,0), cex= 0.8)
    bxp(list(stats= t(sumstats[,c('lW', 'Q1', 'med', 'Q3', 'rW')]), n= sumstats$count), 
        xlab= '', ylab= 'Phred score', main= '', xaxt= 'n', frame.plot= F)
    axis(side= 4, labels= phred2p(seq(2, 35, length.out= 7), 3), at= seq(2, 35, length.out= 7), las= 1)
    mtext(text= 'Base call p-value', side= 4, line= 3, cex= 0.8)

    base_freq<- t(as.matrix(sumstats[,13:17])/sumstats$count)*100
    b<- barplot(base_freq, col= c('limegreen', 'dodgerblue', 'grey24', 'firebrick', 'white'),
        xlab= 'Read position', ylab= 'Base frequency (%)', main= '')
    legend(xpd= T, x= 44, y= 100, legend= c('A', 'C', 'G', 'T', 'N'), 
        fill= c('limegreen', 'dodgerblue', 'grey24', 'firebrick', 'orange'), bty= 'n')

    plot(x= 1:36, base_freq[1,], col= 'limegreen', type='o', pch= 19, lwd= 2, main= '',
        xlab= '', ylab= 'Base frequency (%)', xaxt= 'n', frame.plot= F)
    points(x= 1:36, base_freq[2,], col= 'dodgerblue', type='o', pch= 19, lwd= 2)
    points(x= 1:36, base_freq[3,], col= 'grey24', type='o', pch= 19, lwd= 2)
    points(x= 1:36, base_freq[4,], col= 'firebrick', type='o', pch= 19, lwd= 2)
    points(x= 1:36, base_freq[5,], col= 'orange', type='o', pch= 19, lwd= 1)
    legend(xpd= T, x= 36.7, y= 100, legend= c('A', 'C', 'G', 'T', 'N'), 
        fill= c('limegreen', 'dodgerblue', 'grey24', 'firebrick', 'orange'), bty= 'n')

    axis(side= 1, labels= 1:36, at= 1:36)
    mtext(text= 'Read position', side= 1, line= 2.5, cex= 0.8)
    title(main= 'Base quality and frequency distribution', outer= T)
savePlot('M:/Documents/LabBook/LabBook_Figures/20101207_CAGE_qc.emf', 'emf')

# --------------[ Alignment ]-------------

conn<- odbcConnect(dsn= 'pgCAGE')
aln<- sqlQuery(conn, "
    select rname, 
       count(*) AS n_reads, 
       (count(*)::numeric/(select count(*) from sam_cage_20101207))*100 AS freq 
    from sam_cage_20101207 group by rname
    order by freq desc;
    ")
aln2<- aln
aln[, 1]<- as.character(aln[, 1])
aln[aln$rname=="gi|555853|gb|U13369.1|HSU13369", 1] <- 'hg_rDNA'
barplot(aln$freq, names.arg= aln$rname)

# ------------------[ Tritume ]-----------------------
barplot((sumstats$N_Count/sumstats$count)*100)

p<- 10^(-2/10)

bp<- boxplot(rnorm(n=100))
bxp(bp[1:2])