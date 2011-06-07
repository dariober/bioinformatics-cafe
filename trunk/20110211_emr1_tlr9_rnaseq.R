library(RODBC)
conn<- odbcConnect(dsn= 'pgVitelleschi')

emr1<- sqlQuery(conn, "
    select * from pileup_rnaseq where rname = '2' and pos between 50378276 -100 and 50421628 + 100;
    ")
head(emr1)

par(cex= 0.85, mfrow= c(2,1), mar=c(1,0,0,0), oma= c(4,4,4,2))
plot(emr1$pos, emr1$n_ctrl, type= 'h', col= 'dodgerblue', xaxt= 'n', lwd= 2, xlab= 'Position', ylab= 'No. reads', frame.plot=FALSE)
plot(emr1$pos, emr1$n_lps, type= 'h', col= 'firebrick4', lwd= 2, xlab= 'Position', ylab= 'No. reads', frame.plot=FALSE)
title(main= 'Pig EMR1 2:50378276-50421628', outer=T)

## 
tlr9<- sqlQuery(conn, "
    select * from pileup_rnaseq where rname = '13' and pos between 29049282-100 and 29052128 + 100;
    ")
head(tlr9)

par(cex= 0.85, mfrow= c(2,1), mar=c(1,4,0,0), oma= c(4,0,4,2))
plot(tlr9$pos, tlr9$n_ctrl, type= 'h', col= 'dodgerblue', xaxt= 'n', lwd= 1, xlab= 'Position', ylab= 'No. reads', frame.plot=FALSE)
plot(tlr9$pos, tlr9$n_lps, type= 'h', col= 'firebrick4', lwd= 2, xlab= 'Position', ylab= 'No. reads', frame.plot=FALSE)
title(main= 'Toll-like receptor 9 Precursor (CD289 antigen) Chromosome 13: 29,049,282-29,052,128', outer=T)
