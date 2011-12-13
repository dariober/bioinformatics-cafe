library(edgeR)
library(RODBC)

conn_db<- odbcConnect(dsn= "pgCAGE", case= "nochange", uid= "dberaldi", pwd= "mypassword")

#---------------------------------------------------------------
#                    Transcript analysis
#---------------------------------------------------------------

## Transcripts assembled by cufflinks/cuffcompare denovo.
## Data from 20100512_cuffcompare.sql

raw.data<- sqlFetch(conn_db, 'tmp_denovo_edger')
odbcClose(conn_db)
# -----------[ Prepare for edgeR ]---------------

raw.data[raw.data$fpkm_ctrl == 0 & raw.data$fpkm_lps == 0,]
raw.data[is.na(raw.data)]<- 0 ## Set to 0 count the transcripts found in one lib but not in the other

d<- raw.data[,2:3]
rownames(d) <- raw.data[, 1]
group<- c('ctrl', 'lps')

d <- DGEList(counts = d, group = group)
d <- estimateCommonDisp(d)

# ------------[ Testing DE ]-------------
de.com <- exactTest(d)
names(de.com)
topTags(de.com)

detags.com <- rownames(topTags(de.com)$table) ## See raw counts for top de tags
d$counts[detags.com, ]

hist(de.com$table$p.value)

n.de<- sum(de.com$table$p.value < 0.01) ## Count features with p<0.01
sum(p.adjust(de.com$table$p.value, method = "BH") < 0.05) ## Adjust for multiple testing

top.com <- topTags(de.com, n = n.de)
sum(top.com$table$logFC > 0) ## upregulated in LPS
sum(top.com$table$logFC < 0)

# -------------------------[ Plot DE ]--------------------------
detags500.com <- rownames(topTags(de.com, n = n.de)$table)
plotSmear(d, de.tags = detags500.com, main = "Differential transcript expression in LPS - CTRL\nde novo assembly")
abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100514_edgeR_denovo.bmp', 'bmp')
# ---------------[Send out DE table to CAGE db]-----------------

out.de<- topTags(de.com, n =nrow(de.com$table))$table
out.de<- cbind(dataset_id= '20100514_CTRLvsLPS_denovo', transcript_id= row.names(out.de), out.de)
sqlSave(conn_db, out.de, tablename= 'edger_toptags', rownames= FALSE, append= TRUE)
odbcClose(conn_db)
out.de[1:10,]







sqlSave(conn_db, out.de, tablename= 'tmp_edger_transcript', rownames= FALSE)
