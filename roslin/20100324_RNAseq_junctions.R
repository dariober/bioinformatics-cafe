library(edgeR)
library(RODBC)

conn_db<- odbcConnect(dsn= "pgCAGE", case= "nochange", uid= "dberaldi", pwd= "mypassword")

flag8<- read.table('C:/Tritume/flag8.txt', sep=';', header= T, comment.char = "")
flag2<- read.table('C:/Tritume/flag2.txt', sep=';', header= T, comment.char = "")

par(mfrow= c(2,1), cex= 0.7, mar=c(3,4,0,2))
  h8<- hist(flag8$pos, breaks= c(seq(0, max(flag8$pos), 10^6), max(flag8$pos)), xlab= '', main='')
  legend("topright", legend= 'One read/pair mapped', bty= 'n', cex= 1.5)
  h2<- hist(flag2$pos, breaks= c(seq(0, max(flag2$pos), 10^6), max(flag2$pos)), xlab= 'Chromosome 1', main='')
  legend("topright", legend= 'Both reads/pair mapped', bty= 'n', cex= 1.5)

par(mfrow= c(4,1), cex= 0.7)
  h_lps<- hist(flag8$pos[flag8$dataset_id=='LPS_20100122'], breaks= c(seq(0, max(flag8$pos), 10^6), max(flag8$pos)))
  h_ctrl<- hist(flag8$pos[flag8$dataset_id=='CTRL_20100122'], breaks= c(seq(0, max(flag8$pos), 10^6), max(flag8$pos)))
  h2_lps<- hist(flag2$pos[flag2$dataset_id=='LPS_20100122'], breaks= c(seq(0, max(flag2$pos), 10^6), max(flag2$pos)))
  h2_ctrl<- hist(flag2$pos[flag2$dataset_id=='CTRL_20100122'], breaks= c(seq(0, max(flag2$pos), 10^6), max(flag2$pos)))

# --------------------------[ Junctions ]-----------------------------

raw.data<- read.table('C:/Tritume/junx_count.txt', header=T, sep='\t')
raw.data[is.na(raw.data)]<- 0
raw.data<- raw.data[,1:3] ## Remove gff column

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

# -------------[ Plot DE ]---------------
detags500.com <- rownames(topTags(de.com, n = 500)$table)
plotSmear(d, de.tags = detags500.com, main = "Differential use of splice-junctions in LPS - CTRL")
abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)

# ---------------[Send out DE table to CAGE db]-----------------
out.de<- topTags(de.com, n = nrow(de.com$table))$table
out.de<- cbind(junction_id= row.names(out.de), out.de)
sqlSave(conn_db, out.de, tablename= 'tmp_edger_junction', rownames= FALSE)

#---------------------------------------------------------------
#                    Transcript analysis
#---------------------------------------------------------------

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

raw.data[raw.data$fpkm_ctrl == 0 & raw.data$fpkm_lps == 0,]
raw.data[is.na(raw.data)]<- 0

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
plotSmear(d, de.tags = detags500.com, main = "Differential transcript expression in LPS - CTRL", ylim=c(-10,10))
abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)

# ---------------[Send out DE table to CAGE db]-----------------

out.de<- topTags(de.com, n = nrow(de.com$table))$table
out.de<- cbind(transcript_id= row.names(out.de), out.de)
sqlSave(conn_db, out.de, tablename= 'tmp_edger_transcript', rownames= FALSE)


# --------------------------------[ Tritume ]------------------------------
dim(raw.data)


de.com$table[1:10,]
names(de.com)
names(topTags(de.com, sort.by='logFC'))

topTags(de.com, n= 50)



dim(d)
plot(sort(junx$lps_tpm - junx$ctrl_tpm))

junx[1:10,]
dim(junx)
boxplot(junx[junx$gff_junction==''])



max(flag8$pos)
flag8$pos[flag8$dataset_id=='CTRL_20100122']