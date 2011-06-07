#
# Plot infection rank of the Blackface sheep
#

library(RODBC)
conn<- odbcConnect(dsn= 'pgBlackface')

black<- sqlQuery(conn, 'select * from "Blackface_Sheep"."cBlackF" order by "Rank"')

attach(black)
windows(height= 17/2.54, width= 25/2.54)
par(mfrow=c(3,1), las= 1, mar= c(0, 4, 1, 0), oma= c(4,0,2,0), cex= 0.85)
    b<- barplot(AdTot, col= 'steelblue', ylab= '')
    abline(v= b[10], col= 'blue', lty= 'dotted')
    legend('topleft', legend= 'No. adult worms', cex= 1.2, bty= 'n', bg= 'white')
    barplot(FEC7, col= 'salmon')
    abline(v= b[10], col= 'blue', lty= 'dotted')
    legend('topleft', legend= 'FEC at last time point', cex= 1.2, bty= 'n', bg= 'white')
    barplot(IgA6, col= 'darkolivegreen4')
    abline(v= b[10], col= 'blue', lty= 'dotted')
    legend('topleft', legend= 'IgA at last time point', cex= 1.2, bty= 'n', bg= 'white')
    mtext(text= c(rep(0, 10), (1:47))[seq(1, nrow(black), by= 2)], side= 1, at= as.vector(b)[seq(1, nrow(black), by= 2)], cex= 0.85, xpd= NA)
    mtext(text= 'Infection rank', side= 1, line= 2, xpd= TRUE)
savePlot('M:/Documents/LabBook/LabBook_Figures/20110519_blackface_ranking.emf', 'emf')