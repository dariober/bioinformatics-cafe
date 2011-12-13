library(edgeR)
library(RODBC)

conn_db<- odbcConnect(dsn= "pgCAGE", case= "nochange", uid= "dberaldi", pwd= "mypassword")

#---------------------------------------------------------------
#                    Transcript analysis
#---------------------------------------------------------------

# Transcripts assembled by cufflinks using ensembl GTF file

raw.data<- sqlQuery(conn_db, "
select trans_id.transcript_id, fpkm_ctrl, fpkm_lps FROM
  -- All trascripts found in selected libs, then FPKM for LPS and CTRL libs.
  (select distinct transcript_id from cufflinks_transcript_gtf where source like '20100317_RNAseq_CTRL' or source like '20100317_RNAseq_LPS') as trans_id
  LEFT JOIN
  (select transcript_id, fpkm AS fpkm_ctrl from cufflinks_transcript_gtf where feature like 'transcript' and source like '20100317_RNAseq_CTRL') as ctrl
    ON trans_id.transcript_id = ctrl.transcript_id 
  LEFT JOIN
  (select transcript_id, fpkm AS fpkm_lps from cufflinks_transcript_gtf where feature like 'transcript' and source like '20100317_RNAseq_LPS') as lps
    ON trans_id.transcript_id = lps.transcript_id;
  ")

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
 
n.de<- sum(de.com$table$p.value < 0.01) ## Count features with p<0.01
sum(p.adjust(de.com$table$p.value, method = "BH") < 0.05) ## Adjust for multiple testing

top.com <- topTags(de.com, n = n.de)
sum(top.com$table$logFC > 0) ## upregulated in ctrl
sum(top.com$table$logFC < 0)

# -------------------------[ Plot DE ]--------------------------
detags500.com <- rownames(topTags(de.com, n = 500)$table)
plotSmear(d, de.tags = detags500.com, main = "Differential transcript expression in LPS - CTRL")
abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)

# ---------------[Send out DE table to CAGE db]-----------------

out.de<- topTags(de.com, n = nrow(de.com$table))$table
out.de<- cbind(transcript_id= row.names(out.de), out.de)
sqlSave(conn_db, out.de, tablename= 'tmp_edger_transcript', rownames= FALSE)
