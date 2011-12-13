## Allele specific expression
library(RODBC)
conn_db<- odbcConnect(dsn= "pgCAGE", case= "nochange", uid= "dberaldi", pwd= "mypassword")

snp<- sqlQuery(conn_db,
  "    select 
    lps.rname, lps.pos, lps.allele_1 AS lps_1, lps.count_allele_1 AS lps_count_1, lps.allele_2 AS lps_2, lps.count_allele_2 AS lps_count_2,
    ctrl.allele_1 AS ctrl_1, ctrl.count_allele_1 AS ctrl_count_1, ctrl.allele_2 AS ctrl_2, ctrl.count_allele_2 AS ctrl_count_2,
    log(2, lps.count_allele_1::numeric/lps.count_allele_2) AS lps_logfold, 
    log(2, ctrl.count_allele_1::numeric/ctrl.count_allele_2) AS ctrl_logfold 
  from (select * from tmp_pileup_snp where dataset_id like '20100406_RNAseq_LPS') as lps inner join 
      (select * from tmp_pileup_snp where dataset_id like '20100406_RNAseq_CTRL') as ctrl on
      lps.rname= ctrl.rname and lps.pos= ctrl.pos
  where (lps.allele_1 like ctrl.allele_1) or (lps.allele_2 like ctrl.allele_2)
  order by lps.rname, lps.pos
  ")
par(mfrow= c(1,2), cex= 0.75)
hist(snp$lps_logfold, main= 'Differential allelic expression in\nLPS (red) and CTRL (blue)', xlab= 'Log fold change Allele 1 / Allele 2', border= 'red')
hist(snp$ctrl_logfold, add=T, border= 'blue')
plot(snp$lps_logfold, snp$ctrl_logfold, main='LPS vs CTRL fold change', xlab= 'LPS', ylab= 'CTRL')
abline(a= 0, b= 1, col= 'grey', lwd= 2)
mtext(side= 3, line= -2, text= paste('n=', length(snp$lps_logfold)),adj= 0.05)
names(snp)
dim(snp)