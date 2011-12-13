library(RODBC)
library(plotrix)
conn<- odbcConnect(dsn='pgVitelleschi')
tss_hg_ss_gtf<- sqlFetch(conn, 'fantom5.mphage_mcyte_hg19_on_ss9_ctss_gtf')
tss_ss_gtf<- sqlFetch(conn, 'cage.cage050810_gtf')
tss_hg_ss_gtf[1:10,]
tss_ss_gtf[1:10,]

par(mfrow=c(1,2), cex.axis=0.85, cex.main= 0.95, mar=c(5,2,3,1), oma=c(0,2,0,0))
  plot(density(tss_ss_gtf$tss_exon_dist), type= 'l', col= 'dodgerblue', lwd= 1, frame.plot=F, main= 'Pig CAGE on pig', xlab= '')
  h1<- hist(tss_ss_gtf$tss_exon_dist, freq= FALSE, add=)
  mtext(text= "Distance from exon 5' end\n(+/-1000 bp)", side= 1, line= 3.25, cex= 0.95)
  mtext(text= 'Density', side=2, line=2.5)
  legend('topleft', bty='n', cex= 0.90,
       legend= paste('n tss= ', length(tss_ss_gtf$tss_id), '\nn unique= ', length(unique(tss_ss_gtf$tss_id)), sep=''))
  plot(density(tss_hg_ss_gtf$tss_exon_dist), ylim=c(0, par('usr')[4]), type= 'l', col= 'dodgerblue', lwd= 1, frame.plot=F, main= 'Human CAGE on pig', xlab= '')
  hist(tss_hg_ss_gtf$tss_exon_dist, add=T, freq= T)
  mtext(text= "Distance from exon 5'end\n(+/-1000 bp)", side= 1, line= 3.25, cex= 0.95)
  legend('topleft', bty='n', cex= 0.90,
       legend= paste('n tss= ', length(tss_hg_ss_gtf$tss_id), '\nn unique= ', length(unique(tss_hg_ss_gtf$tss_id)), sep=''))
savePlot('M:/Documents/LabBook/LabBook_Figures/20110110_tss_distance_pig_human.emf', 'emf')

wh1<- weighted.hist(tss_ss_gtf$tss_exon_dist, tss_ss_gtf$tss_tpm, plot=F, breaks= seq(-1000,1000, by= 100))
class(wh1)<- 'histogram'

wh2<- weighted.hist(tss_hg_ss_gtf$tss_exon_dist, tss_hg_ss_gtf$human_tag_count, plot=F, breaks= seq(-1000,1000, by= 100))
class(wh2)<- 'histogram'

par(mfrow=c(1,2), cex.axis=0.85, cex.main= 0.95, mar=c(5,2,3,1), oma=c(0,2,0,0))
  plot(wh1, freq=T, main= 'Pig CAGE on pig\n(weighted)', xlab= '')
  mtext(text= "Distance from exon 5' end\n(+/-1000 bp)", side= 1, line= 3.25, cex= 0.95)
  mtext(text= 'Frequency', side=2, line=2.5, cex= 0.9)
  plot(wh2, freq=T, main= 'Human CAGE on pig\n(weighted)', xlab= '')
  mtext(text= "Distance from exon 5' end\n(+/-1000 bp)", side= 1, line= 3.25, cex= 0.95)
savePlot('M:/Documents/LabBook/LabBook_Figures/20110110_tss_distance_pig_human_weighted.emf', 'emf')

" -------------------[ Transcription level ]------------------- "

tss_hg_ss<- sqlFetch(conn, 'fantom5.mphage_mcyte_hg19_on_ss9_ctss')
tss_ss<- sqlFetch(conn, 'cage.cage050810_tss')
tss_hg_ss[1:10,]

## ACTB
actb_ss<- with(tss_ss, tss_ss[rname == '3' & tss_pos > (3868537 - 100) & tss_pos < 3868537 + 100, c("tss_pos", "tss_tpm")])
actb_hg<- with(tss_hg_ss, tss_hg_ss[rname == '3' & ctss_pos > (3868537 - 100) & ctss_pos < 3868537 + 100, c('ctss_pos', 'human_tpm')])

plot(x= actb_ss$tss_pos, actb_ss$tss_tpm, type= 'h', col= 'dodgerblue', lwd= 2)
points(x= actb_hg$ctss_pos, actb_hg$human_tpm, type= 'h', col= 'firebrick', lwd= 2)

## HPRT
hprt_ss<- with(tss_ss, tss_ss[rname == 'X' & tss_pos > (108419935 - 100) & tss_pos < 108419935 + 100, c("tss_pos", "tss_tpm", 'ctss_id')])
hprt_hg<- with(tss_hg_ss, tss_hg_ss[rname == 'X' & ctss_pos > (108419935 - 100) & ctss_pos < 108419935 + 100, c("ctss_pos", "human_tpm")])

plot(x= hprt_ss$tss_pos, hprt_ss$tss_tpm, type= 'h', col= 'dodgerblue', lwd= 2)
points(x= hprt_hg$ctss_pos, hprt_hg$human_tpm, type= 'h', col= 'firebrick', lwd= 2)

## CSF1R Chromosome 2: 136,526,706-136,552,808 
csf_ss<- with(tss_ss, tss_ss[rname == '2' & tss_pos > 136526706 & tss_pos < 136552808, c("tss_pos", "tss_tpm", 'ctss_id')])
csf_hgss<- with(tss_hg_ss, tss_hg_ss[rname == '2' & ctss_pos > 136526706 & ctss_pos < 136552808, c("ctss_pos", "human_tpm")])

plot(x= csf_ss$tss_pos, csf_ss$tss_tpm, type= 'h', col= 'dodgerblue', lwd= 2, ylim= c(0,10))
points(x= csf_hgss$ctss_pos, csf_hgss$human_tpm, type= 'h', col= 'firebrick', lwd= 2)

## CSF1R in the human sample
csf_hg<- sqlQuery(conn, "select * from macrophage_monocyte_derived_hg19_ctss_bed where rname like 'chr5' and f_start > 149432854 and f_end < 149492935 order by f_start;")
plot(csf_hg$f_start, csf_hg$tpm, type= 'h', col= 'firebrick', lwd= 2)

" ---------------------------------------------------------------- "

mapq<- sqlQuery(conn, 'select mapq from bed_mphage_mcyte_ctss_sscrofa9;')[,1]
mapq[1:10]
-round(quantile(prob= c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 1), -10^(-(mapq)/10)), 4)
quantile(prob= c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 1), mapq)

boxplot(10^(-mapq/10))


