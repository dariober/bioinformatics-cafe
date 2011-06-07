library(RODBC)
library(boot)
conn_db<- odbcConnect(dsn= 'pgCAGE', case='nochange', uid= 'dberaldi', pwd= 'mypassword')
## odbcClose(conn_db)


snp<- sqlFetch(conn_db, 'tmp_snp');
names(snp)
snp[1:10,]

## Test for differential allelic expression by Fisher's test of a contingency table like:
# ------+----+-----
#       | A1 | A2
# ------+----+-----
#  CTRL |  7 | 11
#  LPS  |  3 | 11
# ------+----+-----
#

dae<- vector(length= 0, mode= 'numeric')
for(i in 1:nrow(snp)){
    m<- matrix(c(snp$ctrl_a1[i], snp$ctrl_a2[i], snp$lps_a1[i], snp$lps_a2[i]), byrow= T, nrow= 2)
    rownames(m)<- c('CTRL', 'LPS')
    colnames(m)<- c('A1', 'A2')
    tpval<- fisher.test(m)$p.value ## alternative tests: chisq.test(); prop.test()
    dae<- append(dae, tpval)
    }

## Comapare Fisher's, chisq, prop tests
# dae_prop<- dae
# dae_fish<- dae
# dae_chi<- dae

# length(dae_prop[dae_prop<0.01])
# length(dae_fish[dae_fish<0.01])
# length(dae_chi[dae_chi<0.01])
# pairs(cbind(prop.test= dae_prop, fisher.test= dae_fish, chisq.test= dae_chi), main= 'Comparing tests for detecting differences in proportions')
# savePlot('M:/Documents/LabBook/LabBook_Figures/20100504_RNAseq_SNP_Fig_2.bmp', 'bmp')


## Linear regression Frequency allel 1 CTRL ~ LPS:
## 
names(snp)
## snp_sort<- snp[order(snp$ctrl_af),]
x<- logit(snp$ctrl_af)
y<- logit(snp$lps_af)
model<- lm(x ~ y)
summary(model)
new<- data.frame(x= logit(seq(0.001,0.999, length.out=100)))
ci.model<- data.frame(predict(lm(y~x), new, interval= 'prediction'))
ci.model[1:10,]

## Plot row allele frequency (first allele in alphabetical order) and hist of Fisher's pvalues
par(mfrow=c(1,2), cex= 0.8)
    plot(x= snp$ctrl_af, y= snp$lps_af, xlab= 'CTRL - allele 1', ylab= 'LPS - allele 1', 
        main= 'Allelic expression\nCTRL vs LPS', col= 'darkgrey', pch= 21, cex= 0.5, xlim= c(0,1), ylim= c(0,1))
    points(y= seq(0.01,0.99,length.out= 100), x=inv.logit(sort(ci.model$fit)), type='l', col= 'red', lty= 'dashed', lwd= 2)
    points(y= seq(0.01,0.99, length.out= 100), x= inv.logit(sort(ci.model$lwr)), type='l', col= 'blue', lty= 'dashed', lwd= 2)
    points(y= seq(0.01,0.99, length.out= 100), x= inv.logit(sort(ci.model$upr)), type='l', col= 'blue', lty= 'dashed', lwd= 2)

    hist(dae, main= "Histogram of p-values\n(Fisher's test)", xlab= 'P-value')
savePlot('M:/Documents/LabBook/LabBook_Figures/20100504_RNAseq_SNP_Fig_1.bmp', 'bmp')

## Supposedly differentially expressed alleles
dae_snp<- cbind(snp, fisher_pval= dae)
dae_snp<- dae_snp[order(dae_snp$fisher_pval), ]
write.table(x= dae_snp[dae_snp$fisher_pval<0.01,], file= 'C:/Tritume/snp.txt', row.name= F, sep='\t')

sqlSave(conn_db, dae_snp[dae_snp$fisher_pval<0.01, ], tablename= 'snp_dae', rownames=FALSE)
sqlQuery(conn_db, "COMMENT ON TABLE snp_dae IS 'SNP loci with putative differential allelic expression between RNAseq_CTRL and RNAseq_LPS. Table cretaed by 20100504_RNAseq_SNP.R'")

# ----------------------------[ Tritume ]---------------------------

x<- seq(0,1, length.out=100)
y<- x + rnorm(n= length(x), sd= 0.05)
y[y<0]<- 0
y[y>1]<- 1
plot(x,y)
model<- glm(y~x, family= 'binomial')
ci.model<- data.frame(predict(model, interval= 'prediction', se.fit = TRUE))
ci.model[1:10,]
abline(model)
points(seq(0,1, length.out= nrow(ci.model)), y= ci.model$fit.lwr[order(ci.model$fit.fit)], type= 'l')
points(seq(0,1, length.out= nrow(ci.model)), y= ci.model$fit.upr[order(ci.model$fit.fit)], type= 'l')

ci.lines(model)


m<- matrix(c(385, 103, 8512, 707), nrow=2, byrow=T)
rownames(m)<- c('CTRL', 'LPS')
colnames(m)<- c('A1', 'A2')

n<- matrix(c(800, 850, 8500, 800), nrow=2, byrow=T)

f1<- fisher.test(m)
c1<- chisq.test(m)
names(c1)

Job <- matrix(c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), 4, 4,
dimnames = list(income=c("< 15k", "15-25k", "25-40k", "> 40k"),
                satisfaction=c("VeryD", "LittleD", "ModerateS", "VeryS")))
fisher.test(Job)


TeaTasting <-
matrix(c(3, 1, 1, 3),
       nrow = 2,
       dimnames = list(Guess = c("Milk", "Tea"),
                       Truth = c("Milk", "Tea")))
fisher.test(TeaTasting, alternative = "greater")
## => p=0.2429, association could not be established

# -----------[ Prepare for edgeR ]---------------
raw.data<- data.frame(snp_id= paste(snp$rname, '_', snp$pos, sep=''), snp[, c("ctrl_a1", "lps_a1")])
raw.data[1:10,]

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
plotSmear(d, de.tags = detags500.com, main = "Differential allelic expression in LPS - CTRL")
abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)