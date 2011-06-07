#
# Pig microarray probes corresponding to human genes induced by LPS have been selected 
# to be compared with the corresponding transcripts detected by RNAseq
# 
# See also 20100625_pig_human_affy_rnaseq.sql
#

library(RODBC)
conn_db<- odbcConnect(dsn= 'pgCAGE', case='nochange')
affy<- sqlFetch(conn_db, 'limma_toptable')

m<- matrix(c(24123, 480, 6194, 195), ncol=2)
chisq.test(m)

mr<- matrix(c(9946, 154, 374, 14), ncol=2)
chisq.test(mr)

pig_human<- read.table('C:/Tritume/pig_human_array_affy_rnaseq.txt', header=T, sep='\t')
pig_human[1:10,]
dim(pig_human)

hist(pig_human$fdr_array)
length(unique(pig_human$fdr_array)[unique(pig_human$fdr_array<0.01)])
hist(pig_human$fdr_rnaseq)
length(unique(pig_human$fdr_rnaseq)[unique(pig_human$fdr_rnaseq<0.01)])


hist(affy$adj_p_val)

cor.test(pig_human$logfc_rnaseq[abs(pig_human$logfc_rnaseq)<10], pig_human$logfc_array[abs(pig_human$logfc_rnaseq)<10])
lm1<- lm(pig_human$logfc_rnaseq[abs(pig_human$logfc_rnaseq)<10] ~ pig_human$logfc_array[abs(pig_human$logfc_rnaseq)<10])
summary(lm1)
plot(lm1)

par(mfrow=c(1,2), cex= 0.85)
  boxplot(pig_human[abs(pig_human$logfc_rnaseq)<10, c("logfc_array", "logfc_rnaseq")], 
          names=c("Affymetrix", "RNAseq"), 
          ylab= "LogFC",
          main= 'Spreading of fold change')
  abline(h=c(1,-1), col= 'darkgreen', lty='dotted')
  plot(x= pig_human$logfc_array[abs(pig_human$logfc_rnaseq)<10], 
       y= pig_human$logfc_rnaseq[abs(pig_human$logfc_rnaseq)<10], 
       xlim=c(-2,4.5), 
       ylim=c(-2,4.5),
       xlab= 'Affymetrix LogFC',
       ylab= 'RNAseq LogFC',
       main= 'Correlation in LogFC')
  abline(h=c(1,-1), col= 'darkgreen', lty='dotted')
  abline(v=c(1,-1), col= 'darkgreen', lty='dotted')
  abline(lm1, col='dodgerblue', lwd= 2)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100625_pig_human_affy_rnaseq.emf', 'emf')
