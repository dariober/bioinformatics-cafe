ng<- c(343, 296, 366, 404, 471, 439, 426, 475)
tp<- rep(c(3,5,10,15), each= 2)
tr<- rep(c('LPS', 'CTRL'), 4)

gel<- round(ng*20/1000, 1)

par(cex= 1, mar= c(5,4,4,5))
b1<- barplot(gel, col= rev(rep(heat.colors(4), each= 2)), 
  las= 2, ylab= 'Amount recovered (ug)', ylim= c(0,10),
  main= 'Amount recovered from phen/chlor extraction', cex.main= 0.9)
axis(side= 4, at= seq(0, 500, by= 100)/(1/20*1000), 
  labels= seq(0, 500, by= 100), las= 1, cex= 0.9)
mtext(side=4, line= 3, text= 'Concentration (ng/ul)', cex= 0.9)
mtext(side=1, at=b1, text= tp, line=0.5, cex= 0.85)
mtext(side=1, text= 'Shearing time (min)', line=1.75, cex= 0.9)

lm1<- lm(gel ~ tp)
summary(lm1)
abline(lm1, col= 'dodgerblue', lty= 'dashed', lwd= 2)

savePlot('M:/Documents/LabBook/LabBook_Figures/20101110_chr_shearing.emf', 'emf')
cor.test(gel, tp)


