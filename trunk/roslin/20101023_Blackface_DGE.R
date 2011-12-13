tag_frequency<- read.table('F:/data/20101023_blackface_DGE/tag_frequency.txt', sep='\t', header=T)
tag_frequency[1:10,]
multi_tags<- aggregate(tag_frequency$no_distinct_reads, 
    by= list(tag_frequency$no_distinct_reads), length)

bp<- barplot(multi_tags$x/sum(multi_tags$x), 
             ylab= 'Proportion', xlab= 'Number of different reads',
             cex.axis= 0.9, cex.main= 0.9,
             main= 'Frequency with which different reads map to the same coordinates')

mtext(text= multi_tags$Group.1[c(1, seq(5, length(bp), by= 5))], 
      at= bp[c(1, seq(0, length(bp), by= 5))],
      side= 1, line= 0.5, cex= 0.9
      )
savePlot('F:/data/20101023_blackface_DGE/no_reads_locus.emf', 'emf')
hist(tag_frequency$no_distinct_reads);

est_positions<- read.table('F:/data/20101023_blackface_DGE/est_positions.txt', sep='\t', header=T)
est_positions[1:10,]

hist(est_positions$no_positions, xlab='Number of positions', 
                                 ylab= 'Number of ESTs',
     main='Number of different tags (positions) per EST')
savePlot('F:/data/20101023_blackface_DGE/no_positions_per_EST.emf', 'emf')


# ---------------------[ Difference in library size ]--------------------

lib_size<- read.table('F:/data/20101023_blackface_DGE/counts_per_lib.txt', sep='\t', header=T)

lib_glm<- glm(log(sum) ~ TreatmentGroup, data=lib_size)
summary(lib_glm)

lib_aov<- aov(lib_size$sum ~ lib_size$TreatmentGroup)
summary.lm(lib_aov)
plot(TukeyHSD(lib_aov))
savePlot('F:/data/20101023_blackface_DGE/TukeyHSD.emf', 'emf')
hist(lib_size$sum)


sink(file= 'F:/data/20101023_blackface_DGE/TukeyHSD.txt')
TukeyHSD(lib_aov)
sink()

# ----------------------------[ edgeR differential expression ]-----------

library(edgeR)
library(RODBC)
conn<- odbcConnect('pgBlackface', case= 'nochange')
raw.data<- sqlQuery(conn, 'SELECT * FROM "Blackface_Solexa".out_crosstab_cumcount')
raw.data[is.na(raw.data)]<- 0
rownames(raw.data)<- raw.data[,1]
d<- raw.data[,2:ncol(raw.data)]
d[1:10,]

group<- sqlQuery(conn, 'SELECT * FROM "qrySolexaLambs" ORDER BY "LibraryID"')
group<- group$TreatmentGroup
colnames(raw.data)
group$LibraryID

d <- DGEList(counts = d, group = group)
d <- estimateCommonDisp(d)

de.com <- ?exactTest(d)
names(de.com)
topTags(de.com)

