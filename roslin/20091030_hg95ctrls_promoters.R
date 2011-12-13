library(RODBC)

con<- odbcConnect("pgCAGE", case="nochange")

hg_prom<- sqlFetch(con, "qry_hg95ctrls_promoters")

par(mfrow= c(1,2), cex= 0.85, mar= c(4.5, 4, 1, 1), oma= c(0,0,3,0))
hist(hg_prom$tss_distance, main= "", xlab= "Distance (bp)")
legend("topleft", legend= "5'", bty= "n", cex= 1.5)
legend("topright", legend= "3'      ", bty= "n", cex= 1.5)
hist(hg_prom$tss_distance[hg_prom$tss_distance> -100 & hg_prom$tss_distance< 100],
  main= "", xlab= "Distance (bp)")
title("Distance between CAGE tags (TSS) and 5' ends of neraby genes", outer= TRUE)

# -------------------------------[ Tritume ]-----------------------------------

hg_prom[1:10,]

library("BSgenome")
library(BSgenome.Hsapiens.UCSC.hg19)



getSeq(Hsapiens, names= hg19c$chromosome,   start= hg19c$position - 1,
   end= hg19c$position + 27, strand= hg19c$strand)

