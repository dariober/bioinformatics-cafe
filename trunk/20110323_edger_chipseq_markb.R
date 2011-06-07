#
# Analysis of peak heigth in chip-seq libraries
# See labbook 15/03/2011
# See also 20110321_chipseq_markb.sql
#

source(file.path(Sys.getenv("HOME"), '.Rprofile'))
library(RODBC)
library(edgeR)

conn<- odbcConnect(dsn= 'pgVitelleschi')
## These tables produced by the pipeline findPeaks > mergePeaks > coverageBed.
con_peaks<- sqlQuery(conn, "select rname || '_' || f_start || '_' || f_end || '_' || strand AS peak_id, * from chipseq.peakcoverage where source like 'con2' order by rname, f_start, f_end")
lps_peaks<- sqlQuery(conn, "select rname || '_' || f_start || '_' || f_end || '_' || strand AS peak_id, * from chipseq.peakcoverage where source like 'lps2' order by rname, f_start, f_end")
odbcClose(conn)

con_peaks[1:10,]

# ------------------------------[ Prepare for edgeR ]---------------------------

raw.data<- as.data.frame(cbind(ctrl= con_peaks$read_count, lps= lps_peaks$read_count))
rownames(raw.data)<- con_peaks$peak_id
raw.data[1:10,]

#                        ctrl  lps
# chr1_3395081_3396335_+  130   71
# chr1_4774332_4777028_+ 2037 2103
# chr1_4797077_4799629_+ 1676 1553
# chr1_4846055_4850645_+ 2425 3084
# chr1_5008334_5011146_+  847 1188
# chr1_5072857_5075547_+ 1251 1478

group<- c('ctrl', 'lps')

d <- DGEList(counts = raw.data, group = group)
d <- estimateCommonDisp(d)

# ----------------------------[ Testing DE ]------------------------------------

de.com <- exactTest(d)
names(de.com)
topTags(de.com)

## See raw counts for the top d.e. peaks
detags.com <- rownames(topTags(de.com)$table) 
d$counts[detags.com, ]
 
## Count features with p<0.01
n.de<- sum(de.com$table$p.value < 0.01) 
sum(p.adjust(de.com$table$p.value, method = "BH") < 0.01) ## Adjust for multiple testing

top.com <- topTags(de.com, n = n.de)

names(top.com$table)
sum(top.com$table$logFC > 0) ## upregulated in LPS
sum(top.com$table$logFC < 0)


# -------------------------[ Plot DE ]--------------------------
detags500.com <- rownames(topTags(de.com, n = n.de)$table)
par(mfrow= c(1,1), cex= 0.85)
plotSmear(d, de.tags = detags500.com, main = "Differential precipitation of H3K4me3 in LPS - CTRL")
abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)
savePlot('M:/Documents/LabBook/LabBook_Figures/20110323_chipseq_markb_edger.bmp', 'bmp')

# ---------------[Send out DE table to postgres]-------------------------------

out.de<- topTags(de.com, n = nrow(de.com$table))$table
out.de<- cbind(peak_id= row.names(out.de), out.de)
## Table tmp_edger_transcript copied to 'edger_toptags' with added column 'dataset_id'
## dataset_id for this job is '20100408_LPSvsCTRL_GTF'
sqlSave(conn, out.de, tablename= 'tmp_edger_chipseq', rownames= FALSE)

#                              logConc    logFC PValue FDR
# chr2_32239524_32245164_+   -16.00220 5.597044      0   0
# chr1_6431836_6436135_+     -16.34111 4.633311      0   0
# chr15_31294041_31298563_+  -14.38875 2.892385      0   0
# chr2_34808450_34815090_+   -14.29432 2.830810      0   0
# chr16_84701644_84709262_+  -13.64272 2.176669      0   0
# chr4_144828675_144838303_+ -12.92910 1.813612      0   0



# --------------------------[ Plot individual peaks ]--------------------------

source('U:/Documents/ScriptArchive/R/plot_template_single.R')
setwd('F:/Tritume')

plot_peaks<- function(rname, from, to, peaks){
    pp<- import.pileup(c('F:/Tritume/con2_k4.sorted.bam', 'F:/Tritume/lps2_k4.sorted.bam'), rname=rname, from= from-10000, to= to+10000, debug= TRUE)
    plot(y= pp$count.1[seq(1, nrow(pp), by= 100)], x= pp$pos[seq(1, nrow(pp), by= 100)], type='l', ylim=c(0, max(pp[,4:5])), col= 'blue', lwd= 2, xlab= 'Position', ylab= 'Count', main= peak )
    points(y= pp$count.2[seq(50, nrow(pp), by= 100)], x= pp$pos[seq(50, nrow(pp), by= 100)], type='l', col= 'red', lwd= 2)
    abline(v= c(from, to), lty= 'dotted', col= 'grey40', lwd= 2)
    legend('topright', legend=c('ctrl', 'lps'), col= c('blue', 'red'), lwd= 2)
    }
rname<- 'chr2'
from<- 32239524
to<- 32245164
peak<- 'chr2_32239524_32245164_+'
plot_peaks(rname, from, to, peaks)

rname<- 'chr1'
from<- 6431836
to<- 6436135
peak<- 'chr1_6431836_6436135_+'
plot_peaks(rname, from, to, peaks)

peak<- 'chr3_96361177_96368399_+'
rname<- 'chr3'
from<- 96361177
to<- 96368399
plot_peaks(rname, from, to, peaks)

peak<- 'chr15_31294041_31298563_+'
rname<- 'chr15'
from<- 31294041
to<- 31298563
plot_peaks(rname, from, to, peaks)

peak<- "chr7_105638041_105642799_+"
rname<- "chr7"
from<- 105638041
to<- 105642799
plot_peaks(rname, from, to, peaks)

# -------------------[ DESeq Differential precipitation ]----------------------

library(DESeq)

conn<- odbcConnect(dsn= 'pgVitelleschi')
## These tables produced by the pipeline findPeaks > mergePeaks > coverageBed.
con_peaks<- sqlQuery(conn, "select rname || '_' || f_start || '_' || f_end || '_' || strand AS peak_id, * from chipseq.peakcoverage where source like 'con2' order by rname, f_start, f_end")
lps_peaks<- sqlQuery(conn, "select rname || '_' || f_start || '_' || f_end || '_' || strand AS peak_id, * from chipseq.peakcoverage where source like 'lps2' order by rname, f_start, f_end")
odbcClose(conn)

con_peaks[1:10,]

# ------------------------------[ Prepare for DESeq ]---------------------------
## Same as edger
countsTable<- as.data.frame(cbind(ctrl= con_peaks$read_count, lps= lps_peaks$read_count))
rownames(countsTable)<- con_peaks$peak_id
countsTable[1:10,]

#                        ctrl  lps
# chr1_3395081_3396335_+  130   71
# chr1_4774332_4777028_+ 2037 2103
# chr1_4797077_4799629_+ 1676 1553
# chr1_4846055_4850645_+ 2425 3084
# chr1_5008334_5011146_+  847 1188
# chr1_5072857_5075547_+ 1251 1478

conds<- c('ctrl', 'lps')

# ------------------------------------------------------------------------------
# DESeq Differential Expression
# ------------------------------------------------------------------------------

cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )       ## Factor to adjust library sizes to account for differences in seq depth

## For each gene, use the mean across all samples to estimate the associated variance
## method= 'blind' allows for no replication.
cds <- estimateVarianceFunctions( cds, method= 'blind' )

res <- nbinomTest( cds, "ctrl", "lps")  ## Binomial test for D.E.
head(res)

## Send to postgres
sqlQuery(conn, "set search_path to 'chipseq'")
sqlSave(conn, res, tablename= 'deseq_nbinomtest', rownames= FALSE)
sqlQuery(conn, "COMMENT ON TABLE deseq_nbinomtest IS 'Output of DESeq from comparing ctrl/lps. See labbook 24/03/2011 and 20110323_edger_chipseq_markb.R. Dataset is the same as for edger.' ")

:> head(res)
#                      id  baseMean baseMeanA  baseMeanB foldChange log2FoldChange      pval padj resVarA resVarB
#1 chr1_3395081_3396335_+  102.7475  139.1761   66.31884  0.4765101    -1.06942140 0.1710383    1      NA      NA
#2 chr1_4774332_4777028_+ 2072.5642 2180.7832 1964.34522  0.9007522    -0.15079784 0.5572685    1      NA      NA
#3 chr1_4797077_4799629_+ 1622.4547 1794.3017 1450.60776  0.8084525    -0.30676503 0.2705804    1      NA      NA
#4 chr1_4846055_4850645_+ 2738.4182 2596.1704 2880.66603  1.1095828     0.15001731 0.5486946    1      NA      NA
#5 chr1_5008334_5011146_+ 1008.2295  906.7861 1109.67291  1.2237427     0.29130025 0.4413344    1      NA      NA
#6 chr1_5072857_5075547_+ 1359.9277 1339.3028 1380.55266  1.0307995     0.04376377 0.9058069    1      NA      NA
:> 

plot( res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( res$padj < .1, "red", "black" ),
     ylab= 'Log2 fold change', xlab= 'Base mean', cex.main= 0.85, main= 'Peaks found D.E. (padj <.1) in\nrelation to fold change and base mean'
abline(h= c(-2, 2), col= 'dodgerblue', lwd= 2)
savePlot('M:/Documents/LabBook/LabBook_Figures/20110323_chipseq_markb_deseq.bmp', 'bmp')

resSig <- res[ res$padj < 0.1, ]
dim(resSig)

# -------------------------[ EdgeR vs DESeq ]-----------------------------------

## This is a peak found DE by edgeR (fdr 10e-8) but not by DESeq

peak<- "chr1_4797077_4799629_+"
rname<- "chr1"
from<- 4797077
to<- 4799629
plot_peaks(rname, from, to, peaks)
