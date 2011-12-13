library(RODBC)
conn<- odbcConnect(dsn= 'pgVitelleschi')
odbcClose(conn)
lib_expr<- sqlQuery(conn, "
  select min(treatment_group) AS treatment_group, tag_expression.library_id, sum(tag_count) AS total_expression, count(tag_id) AS no_tags 
  from tag_expression inner join treatment_groups on tag_expression.library_id = treatment_groups.library_id
  group by tag_expression.library_id
  order by treatment_group, total_expression;
  ")

par(mar=c(5,6,2,1))
barplot(lib_expr$total_expression/1000, 
  col= rep(c('lightgrey', 'dodgerblue', 'firebrick'), each= 5),
  names.arg= lib_expr$library_id, las= 2, ylab= '')
mtext(text='Total tag count\nx1000', side=2, line= 4)
legend('topleft', legend= c('Control', 'Resistant', 'Susceptible'), 
  fill= c('lightgrey', 'dodgerblue', 'firebrick'), bty= 'n')
savePlot('M:/Documents/LabBook/LabBook_Figures/20101224_blackface_tag_counts.emf', 'emf')

tpm<- sqlQuery(conn, "
  select tag_expression.tag_id,
     treatment_groups.treatment_group, 
     sum(tag_tpm) AS group_tpm
  from treatment_groups inner join tag_expression on tag_expression.library_id = treatment_groups.library_id
  group by treatment_groups.treatment_group, tag_expression.tag_id
  ")
tpm[1:10,]
hist(log(tpm$group_tpm), xlab='log(tpm)', main='Histogram of tags per million (TPM) per treatment group\n(TPM summed within treatment groups)')
abline(v= log(2), lty='dotted', col='blue', lwd=2)
savePlot('M:/Documents/LabBook/LabBook_Figures/20101224_blackface_hist_tpmBYgroup.emf', 'emf')


library(edgeR)
raw.data<- sqlFetch(conn, 'sig_tags_ct')
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
names(de.com$table)
topTags(de.com.RS, n= 30) ## nrow(raw.data)

## concatenate tables of DE expression:
edger_toptags<- rbind(
      cbind(pair= 'RS', topTags(de.com.RS, n= nrow(de.com.RS$table))$table),
      cbind(pair= 'CR', topTags(de.com.RS, n= nrow(de.com.CR$table))$table),
      cbind(pair= 'CS', topTags(de.com.RS, n= nrow(de.com.CS$table))$table)
      )

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

# ---------------------------[ Tritume ]-------------------------

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