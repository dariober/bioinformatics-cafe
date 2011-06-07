library(RODBC)
library(reshape)

## Import read quality (phred converted) for Test8
conn_db<- odbcConnect(dsn= "pgCAGE", case= "nochange", uid= "dberaldi", pwd= "mypassword")

## Select only the insert of each CAGE construct, not the linker AGACAGCAG
qscore<- sqlQuery(channel= conn_db, query= 
  "SELECT 
     solexa_to_phred(substring(quality from 10 for (length(quality) - 9))) AS phred
   FROM fastq 
   WHERE fastq_name like 'test8'
   ORDER BY random()
   LIMIT 100000")

odbcClose(conn_db)

qscore$phred<- as.character(qscore$phred)

## Make each base a column of a dataframe
qscore<- colsplit(qscore$phred, ",", names= paste("base", seq(1, 27), sep= ""))

## Average of quality for each read
meanq<- apply(qscore, 1, mean)
hist(meanq, xlab= "Phred quality score", main= "Histogram of read average quality")

## Distribution of quality sccore base by base
par(mfrow= c(5,6), mar= c(1,1,1,1), oma= c(2,2,2,2))
for(i in 1:dim(qscore)[2]){
  hist(qscore[, i], xlab= "", main= "",
    ylim= c(0, 35000))
  mtext(paste("Base no.", i), side= 3, line= -2.5, cex= 0.75)
  }

#--------------------------------------[ Tritume ]-----------------------------

plot(as.numeric(qscore[1,]), ylim= c(2, 34), type= "b")  
abline(h= 20, lty= "dotted", col= "blue")
points(as.numeric(qscore[3,]), type= "b")

boxplot(qscore[1:10000,])


meanq[1:10]

qscore[1:10,]

dim(qscore)
object.size(qscore)
memory.size()
memory.limit()


boxp<- seq(1:100)
bp<- boxplot(boxp)
points(rep(1, length(boxp)), boxp)

boxp<- c(rnorm(n= 50), rnorm(n= 50)+10)
bp<- boxplot(boxp)
points(rep(1, length(boxp)), boxp)
hist(boxp)