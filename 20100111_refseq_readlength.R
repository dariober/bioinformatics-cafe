library(RODBC)

conn_db<- odbcConnect("pgCAGE", case= "nochange", "dberaldi", "mypassword")

s7_1<- sqlQuery(conn_db, "select *, length(quality) from tmp_pe_s_7_1 limit 100000 offset 100000")

s8_2<- sqlQuery(conn_db, "select *, length(quality) from tmp_pe_s_8_2 limit 100000 offset 100000")

odbcClose(conn_db)

summary(s7_1$length)
summary(s8_2$length)

par(mfrow= c(1,2), oma= c(0,0,2.5,0), cex= 0.75)
  hist(s7_1$length, 
  main= "s_7_1 (LPS)", 
  xlab= "Read length")
  
  hist(s8_2$length, 
  main= "s_8_2 (CTRL)", 
  xlab= "Read length")
title(main= "Read length after filtering\nSample of 100000 reads",
  outer= TRUE)
  

