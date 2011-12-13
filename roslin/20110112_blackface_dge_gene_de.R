library(RODBC)
library(edgeR)

conn<- odbcConnect(dsn= 'pgVitelleschi')
# odbcClose(conn)

raw.data<- sqlFetch(conn, 'gene_coverage_ct')
treatment<- sqlQuery(conn, 'SELECT * FROM treatment_groups ORDER BY library_id')
raw.data[is.na(raw.data)]<- 0 ## Set to 0 count the transcripts found in one lib in another

d<- raw.data[,2:ncol(raw.data)]
rownames(d) <- raw.data[, 1]
## Make sure names in the data-table are the same as in the target (treatment) table:
all(names(d) == treatment$library_id)

group<- treatment$treatment_group

d <- DGEList(counts = d, group = group)
d <- estimateCommonDisp(d)

# ------------[ Testing DE ]-------------
de.com.RS <- exactTest(d, pair= c('R', 'S'))
de.com.CR <- exactTest(d, pair= c('C', 'R'))
de.com.CS <- exactTest(d, pair= c('C', 'S'))
topTags(de.com.RS, n= 30) ## nrow(raw.data)

## concatenate tables of DE expression:
edger_toptags<- rbind(
      cbind(pair= 'RS', topTags(de.com.RS, n= nrow(de.com.RS$table))$table),
      cbind(pair= 'CR', topTags(de.com.CS, n= nrow(de.com.CS$table))$table),
      cbind(pair= 'CS', topTags(de.com.CR, n= nrow(de.com.CR$table))$table)
      )
## Upload to postrges
sqlSave(conn, edger_toptags, row.names=T)
sqlQuery(conn, 'ALTER TABLE edger_toptags RENAME COLUMN rownames TO gene_id;')
sqlQuery(conn, "comment on table edger_toptags is 'Differential expression analysis done at the level of gene_id. See 20110112_blackface_dge_gene.R. Level of gene expression calculated by HTSeq-count';")

# -----------------------[ Stub: GO analysis ]-------------------

library(RODBC)
library(goseq) ## Make sure you have downloaded formatted a database of GO terms
library(GO.db)

conn<- odbcConnect(dsn='pgVitelleschi')
tested<- sqlQuery(conn, "
  SELECT * FROM edger_toptags WHERE pair like 'RS'
  ")
tested[1:10,]
hist(tested$logconc)
# genes<- as.integer(p.adjust(tested$pvalue[tested$logfc != 0], method = "BH") < 0.05)
genes<- as.integer(abs(tested$logfc) > 1.5 & tested$logconc > -30)

names(genes)<- tested$gene_id
table(genes)

pwf<- nullp(genes, "bosTau4", "ensGene")
GO.wall<- goseq(pwf, "bosTau4", "ensGene")
head(GO.wall)
names(GO.wall)

enriched.GO<- GO.wall$category[GO.wall$upval < 0.05]

for (go in enriched.GO) {
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
  }

genes[1:10]


# ---------------------------[ Tritume ]-------------------------


hist(de.com.RS$table$p.value)

detags.com <- rownames(topTags(de.com)$table) ## See raw counts for top de tags
d$counts[detags.com, ]
 
n.de<- sum(de.com.RS$table$p.value < 0.01) ## Count features with p<0.01
sum(p.adjust(de.com.RS$table$p.value, method = "fdr") < 0.05) ## Adjust for multiple testing

top.com <- topTags(de.com.RS, n = n.de)

names(top.com$table)
sum(top.com$table$logFC > 0) ## upregulated in LPS
sum(top.com$table$logFC < 0)

# -------------------------[ Plot DE ]--------------------------
detags500.com <- rownames(topTags(de.com, n = n.de)$table)
par(mfrow= c(1,2))
    plotSmear(d, de.tags = detags500.com, main = "Differential transcript expression in LPS - CTRL")
    abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)

    ## Zoom to remove low expressed transcripts
    plotSmear(d, de.tags = detags500.com, main = "Differential transcript expression in LPS - CTRL", xlim= c(-25, -5), ylim= c(-7,7))
    abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)
## savePlot('M:/Documents/LabBook/LabBook_Figures/', 'bmp')


raw.data[which(raw.data$tag_id == '20_66278785_0'), ]

rs<- as.matrix(raw.data[,2:ncol(raw.data)])
rsum<- rowSums(rs)

summary(rsum)
raw.data[which(rowSums(raw.data[,2:ncol(raw.data)]) == 4), ]


expr_mat<- raw.data[,2:ncol(raw.data)]
rownames(expr_mat) <- raw.data[,1]
expr_mat[1:10,1:10]

corr_mat<- cor(expr_mat)
heatmap(corr_mat, margins=c(10,10))