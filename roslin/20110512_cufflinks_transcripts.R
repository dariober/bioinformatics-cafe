# ------------------------------------------------------------------------------
# See Labbook 12/05/2011
# ------------------------------------------------------------------------------

library(RODBC)
conn<- odbcConnect('pgVitelleschi')

sqlQuery(conn,
    "select cross_tab($$ select transcript_id,  source, fpkm from cufflinks_transcript_gtf where feature = 'transcript' $$, 'cufflinks_ct')"
)
cuff<- sqlQuery(conn, 'select * from cufflinks_ct')
sqlQuery(conn, "drop table cufflinks_ct")
row.names(cuff)<- cuff[,1]
cuff<- cuff[,-1]
names(cuff)<- c('am_ctrl', 'am_lps', 'bmdm_ctrl', 'bmdm_lps')
cuff<- as.matrix(cuff)
cuff.r<- cuff[rowMeans(cuff[,-1])> 2^5,]
cuff.r[1:10,]
dim(cuff.r)

##
par(cex= 0.85)
hcuff<- hclust(dist(t(cuff.r)))
plot(hcuff)
cor(cuff.r)

## Remove the 10% lowest expressed over the all dataset.
hist(log2(rowMeans(cuff)), breaks= 50)
exprq<- quantile(log2(as.vector(cuff)[as.vector(cuff)>0]), seq(0,1, by= 0.1))
q10<- exprq[2]
cuffq<- cuff[rowSums(cuff)>q10,]


hist(log2(as.vector(cuffq)[as.vector(cuffq)>0]), breaks= 50)

boxplot(log2(cuff.r))
summary(cuff)
length(as.vector(cuff))

2^0
log2(16)