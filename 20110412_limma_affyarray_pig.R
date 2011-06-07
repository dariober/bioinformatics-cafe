# -----------------------------------------------------------------------------
# Microarray analysis of pig BMDM arrays using LIMMA
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Normalization
# -----------------------------------------------------------------------------

library(affy)
celfile_dir<- "U:/Documents/affymetrix/DATA" ## Might need to unzip first.
setwd(celfile_dir)
Data <- ReadAffy() ##read data in working directory
eset <- rma(Data)

# -----------------------------[ Export ]--------------------------------------

## This should be the same as table affymetrix.affyarray_design in postgres
array_design<- read.table('U:/Documents/affymetrix/pig_bmdm_array_design.txt', sep= '\t', header= TRUE)
#:> array_design
#               slide_id  pig time_point cell_culture_replicate                              cel_file
#1   s001_LWL1_0h_repl_1 LWL1          0                      1 1741_388_Kapetonovic_1783_388_001.CEL
#2   s002_LWL1_2h_repl_1 LWL1          2                      1 1742_388_Kapetonovic_1784_388_002.CEL
#3   s003_LWL1_7h_repl_1 LWL1          7                      1 1743_388_Kapetonovic_1785_388_003.CEL
#4  s004_LWL1_24h_repl_1 LWL1         24                      1 1744_388_Kapetonovic_1786_388_004.CEL
#5   s005_LWL1_0h_repl_2 LWL1          0                      2 1745_388_Kapetonovic_1787_388_005.CEL
#6   s006_LWL1_2h_repl_2 LWL1          2                      2 1746_388_Kapetonovic_1788_388_006.CEL
#7   s007_LWL1_7h_repl_2 LWL1          7                      2 1747_388_Kapetonovic_1789_388_007.CEL
#8  s008_LWL1_24h_repl_2 LWL1         24                      2 1748_388_Kapetonovic_1790_388_008.CEL
#9   s009_LWL4_0h_repl_1 LWL4          0                      1 1749_388_Kapetonovic_1791_388_009.CEL
#10  s010_LWL4_2h_repl_1 LWL4          2                      1 1750_388_Kapetonovic_1792_388_010.CEL
#11  s011_LWL4_7h_repl_1 LWL4          7                      1 1751_388_Kapetonovic_1793_388_011.CEL
#12 s012_LWL4_24h_repl_1 LWL4         24                      1 1752_388_Kapetonovic_1794_388_012.CEL
#13  s013_LWL4_0h_repl_2 LWL4          0                      2 1753_388_Kapetonovic_1795_388_013.CEL
#14  s014_LWL4_2h_repl_2 LWL4          2                      2 1754_388_Kapetonovic_1796_388_014.CEL
#15  s015_LWL4_7h_repl_2 LWL4          7                      2 1755_388_Kapetonovic_1797_388_015.CEL
#16 s016_LWL4_24h_repl_2 LWL4         24                      2 1756_388_Kapetonovic_1798_388_016.CEL
#17  s017_LWL5_0h_repl_1 LWL5          0                      1 1757_388_Kapetonovic_1799_388_017.CEL
#18  s018_LWL5_2h_repl_1 LWL5          2                      1 1758_388_Kapetonovic_1813_388_018.CEL
#19  s019_LWL5_7h_repl_1 LWL5          7                      1 1759_388_Kapetonovic_1801_388_019.CEL
#20 s020_LWL5_24h_repl_1 LWL5         24                      1 1760_388_Kapetonovic_1802_388_020.CEL
#21  s021_LWL5_0h_repl_2 LWL5          0                      2 1761_388_Kapetonovic_1803_388_021.CEL
#22  s022_LWL5_2h_repl_2 LWL5          2                      2 1762_388_Kapetonovic_1804_388_022.CEL
#23  s023_LWL5_7h_repl_2 LWL5          7                      2 1763_388_Kapetonovic_1805_388_023.CEL
#24 s024_LWL5_24h_repl_2 LWL5         24                      2 1764_388_Kapetonovic_1814_388_024.CEL


# Convert expression matrix to 'long' format and send to postgres
exprs_pig<- exprs(eset)

## Rename columns to have the slide_id rather than the name of the cel file
slide_id<- as.character(array_design$slide_id[match(colnames(exprs_pig), array_design$cel_file)])
colnames(exprs_pig)<- slide_id

## Write matrix to file
write.table(x= exprs_pig, file= 'U:/Documents/affymetrix/20110412_affymetrix_pig_rma.txt', sep= '\t', quote= FALSE)

## Convert 
exprs_pig_probes<- rep(rownames(exprs_pig), ncol(exprs_pig))
exprs_pig_arrays<- rep(colnames(exprs_pig), each= nrow(exprs_pig))
values<- as.vector(exprs_pig)
pig_array<- data.frame(probe_set_id= exprs_pig_probes, slide_id= exprs_pig_arrays, log2_intensity= values)

library(RODBC)
conn<- odbcConnect(dsn= 'pgVitelleschi')
sqlSave(conn, pig_array, tablename= 'affyarray_bmdm', rownames=FALSE)  
sqlQuery(conn, 'ALTER TABLE affyarray_bmdm SET SCHEMA affymetrix');

library(limma)

# -------------------------[ Import original dataset ]-------------------------

pig_array_ori<- read.table('U:/Documents/affymetrix/20110412_affymetrix_pig_rma.txt', sep='\t', header=TRUE, stringsAsFactors= FALSE)

colnames(pig_array_ori)[1:10]
rownames(pig_array_ori)[1:10]

## Input should look like this:
## --------------------------------------------------------------
# :> pig_array_ori[1:5,1:5]
#                s001_LWL1_0h_repl_1 s002_LWL1_2h_repl_1 s003_LWL1_7h_repl_1 s004_LWL1_24h_repl_1 s005_LWL1_0h_repl_2
# AFFX-BioB-3_at            8.249044            8.501810            8.496941             8.726919            8.400928
# AFFX-BioB-5_at            8.125517            8.366011            8.510154             8.691518            8.220821
# AFFX-BioB-M_at            8.381620            8.782618            8.844169             9.078234            8.695613
# AFFX-BioC-3_at            9.772719            9.957900            9.953191            10.206591            9.928164
# AFFX-BioC-5_at            9.087081            9.282562            9.338154             9.635377            9.253312
## --------------------------------------------------------------

## Average replicates

limma_affy<- cbind(
    LWL1_0h=  rowMeans(pig_array_ori[, c(1, 5)]),      
    LWL1_2h=  rowMeans(pig_array_ori[, c(2, 6)]),
    LWL1_7h=  rowMeans(pig_array_ori[, c(3, 7)]),
    LWL1_24h= rowMeans(pig_array_ori[, c(4, 8)]),
    LWL4_0h=  rowMeans(pig_array_ori[, c(9, 13)]),
    LWL4_2h=  rowMeans(pig_array_ori[, c(10, 14)]),
    LWL4_7h=  rowMeans(pig_array_ori[, c(11, 15)]),
    LWL4_24h= rowMeans(pig_array_ori[, c(12, 16)]),
    LWL5_0h=  rowMeans(pig_array_ori[, c(17, 21)]),
    LWL5_2h=  rowMeans(pig_array_ori[, c(18, 22)]),
    LWL5_7h=  rowMeans(pig_array_ori[, c(19, 23)]),
    LWL5_24h= pig_array_ori[, c(20)]    ## Note: Array in column 25 "s024_LWL5_24h_repl_2" should be excluded (see LabBook 09/07/2010)
    )

limma_affy[1:10,]

# -----------------------------------------------------------------------------
# Cluster arrays
# -----------------------------------------------------------------------------
# Show that "s024_LWL5_24h_repl_2" is not right!
d<- dist(t(pig_array_ori))           ## Produce dissimilarity matrix
clus_d<- hclust(d)    ## Produce hierarchical cluster
par(mar= c(1,1,3,12))
plot(as.dendrogram(clus_d), horiz= TRUE, axes= FALSE, main= 'Clusters from the 24 pig arrays')          ## Plot it
savePlot('M:/Documents/LabBook/LabBook_Figures/20110412_bmdm_pig_arrays_clusters.emf', 'emf')
## This figure re-formatted in PowerPoint.

# -----------------------------------------------------------------------------
# Import or prepare target file
# -----------------------------------------------------------------------------

# NB: It is essential that the arrays in the *target* dataframe are in the same order as in the inetnsity *matrix* <=== 
## targets_ori2<- read.table('U:/Documents/affymetrix/pig_bmdm_array_design_avg.txt', sep='\t',header=T)
pig<- rep(c('LWL1', 'LWL4', 'LWL5'), each= 4)
time_point<- rep(c(0,2,7,24), times= 3)
slide_id<- paste(pig, paste(time_point, 'h', sep= ''), 'repl_1', sep= '_') ## Formatted for consistency with previous versions of this script.
targets_ori<- data.frame(slide_id, pig, time_point= as.integer(time_point))

## Target file should be like this:
#
#           slide_id  pig time_point
# 1   LWL1_0h_repl_1 LWL1          0
# 2   LWL1_2h_repl_1 LWL1          2
# 3   LWL1_7h_repl_1 LWL1          7
# 4  LWL1_24h_repl_1 LWL1         24
# 5   LWL4_0h_repl_1 LWL4          0
# 6   LWL4_2h_repl_1 LWL4          2
# 7   LWL4_7h_repl_1 LWL4          7
# 8  LWL4_24h_repl_1 LWL4         24
# 9   LWL5_0h_repl_1 LWL5          0
# 10  LWL5_2h_repl_1 LWL5          2
# 11  LWL5_7h_repl_1 LWL5          7
# 12 LWL5_24h_repl_1 LWL5         24

# -----------------------------------------------------------------------------
# Prepare design matrix
# -----------------------------------------------------------------------------

design<- model.matrix(~0 + targets_ori$pig)
rownames(design)<- paste(targets_ori$pig, targets_ori$time_point, sep= '.')
design<- cbind(design, tp.2vs0= c(0,1,0,0, 0,1,0,0, 0,1,0,0))
design<- cbind(design, tp.7vs0= c(0,0,1,0, 0,0,1,0, 0,0,1,0))
design<- cbind(design, tp.24vs0= c(0,0,0,1, 0,0,0,1, 0,0,0,1))

## Note: The same matrix is produced by (not tested):
## design<- model.matrix(~0 + targets_ori$pig + as.factor(targets_ori$time_point))
## 

# :> design
#         pig_idLWL1 pig_idLWL4 pig_idLWL5 tp.2vs0 tp.7vs0 tp.24vs0
# LWL1.0           1          0          0       0       0        0
# LWL1.2           1          0          0       1       0        0
# LWL1.7           1          0          0       0       1        0
# LWL1.24          1          0          0       0       0        1
# LWL4.0           0          1          0       0       0        0
# LWL4.2           0          1          0       1       0        0
# LWL4.7           0          1          0       0       1        0
# LWL4.24          0          1          0       0       0        1
# LWL5.0           0          0          1       0       0        0
# LWL5.2           0          0          1       1       0        0
# LWL5.7           0          0          1       0       1        0
# LWL5.24          0          0          1       0       0        1

## Note: Columns 4,5,6 will contrast 0 vs 2, 7, 24, respectively.
## Why? Because the rows with 0 (LWL1.0, LWL4.0, LWL4.0) do not have 1.
## If you remove columns 5 and 6, time point 2 will be compared against all
## the other time points (logfc will be [logfc= 2 - avg(0,7,24)] ).

# -----------------------------------------------------------------------------
# Differential expression
# -----------------------------------------------------------------------------

fit <- lmFit(limma_affy, design)
fit <- eBayes(fit)
topTable(fit, coef= 4)

limma_toptable<- topTable(fit, coef= 4)[0,]
limma_toptable$coef<- character()
coefs<- c('2vs0', '7vs0', '24vs0')
for (i in 1:3){
    tt<- topTable(fit, coef= i+3, number= nrow(limma_affy))
    tt$coef<- coefs[i]
    limma_toptable<- rbind(limma_toptable, tt)
}

# -----------------------------------------------------------------------------
# Export diff. exprs. and summary plots
# -----------------------------------------------------------------------------

sqlSave(conn, limma_toptable, tablename= 'limma_toptable', rownames= FALSE)
sqlQuery(conn, 'ALTER TABLE limma_toptable SET SCHEMA affymetrix');
sqlQuery(conn, "COMMENT ON TABLE affymetrix.limma_toptable IS 'Output of LIMMA function topTable for pig arrays BMDM. See 20110412_limma_affyarray_pig.R, labbook 14/04/2011.'");

comparison<- '2vs0'
comparison<- '7vs0'
comparison<- '24vs0'
plotname<- file.path('M:/Documents/LabBook/LabBook_Figures', paste('20110412_diff_expr_affymetrix', comparison, 'tiff', sep= '.'))
tiff(plotname,
      res= 150,                             ## Resolution in dpi
      pointsize= 8,                         ## This determines the size of the writings
      units= 'cm', width = 10, height = 5.5 ## Unit of measure and size
      )
par(mfrow=c(1,2), cex=0.8, mar= c(4,3,0,1), oma= c(0,0,4,0), mgp= c(2.1,0.8,0), fg= 'grey25', col.axis= 'grey25')
hist(limma_toptable$P.Value[limma_toptable$coef == comparison], main= '', xlab='p-value')
par(las= 1)
plot(limma_toptable$AveExpr[limma_toptable$coef == comparison],
    limma_toptable$logFC[limma_toptable$coef == comparison],
    col= ifelse(limma_toptable$adj.P.Val[limma_toptable$coef == comparison]<0.01, 'blue', 'grey25'),
    pch=16, cex=0.2, 
    main= '', xlab= 'Log average expression', ylab= 'log2 fold change')
title(outer= T, main= paste("Expression differences at", comparison) )
dev.off()


# -----------------------------------------------------------------------------
# Venn diagrams
# -----------------------------------------------------------------------------

#
# Sort table to have probes in different comparisons in the same order:
#
limma_toptable<- limma_toptable[order(limma_toptable$coef, limma_toptable$ID), ]
limma_toptable[1:100,]

## Up-regulated
degenes<- limma_toptable[limma_toptable$coef == '2vs0',]
A<- ifelse(degenes$logFC>0 & degenes$adj.P.Val<0.01, 1, 0)
degenes<- limma_toptable[limma_toptable$coef == '7vs0',]
B<- ifelse(degenes$logFC>0 & degenes$adj.P.Val<0.01, 1, 0)
degenes<- limma_toptable[limma_toptable$coef == '24vs0',]
C<- ifelse(degenes$logFC>0 & degenes$adj.P.Val<0.01, 1, 0)

D<- list()
D$table<- table(A, B, C)
D$labels<- c("2 vs 0","7 vs 0","24 vs 0")

## Show counts in each set:
vC<- vennCounts(cbind('2vs0'= A, '7vs0'= B, '24vs0'= C))
vC_perc<- vC[,4] / sum(vC[,4]) * 100
vC<- cbind(as.matrix(vC), vC_perc)

#> vC
#     2vs0 7vs0 24vs0 Counts    vC_perc
#[1,]    0    0     0  20780 86.1418563
#[2,]    0    0     1    860  3.5650624
#[3,]    0    1     0   1002  4.1537122
#[4,]    0    1     1    881  3.6521162
#[5,]    1    0     0    108  0.4477055
#[6,]    1    0     1     29  0.1202172
#[7,]    1    1     0    143  0.5927953
#[8,]    1    1     1    320  1.3265348
#> 


## Down-regulated
degenes<- limma_toptable[limma_toptable$coef == '2vs0',]
a<- ifelse(degenes$logFC<0 & degenes$adj.P.Val<0.01, 1, 0)
degenes<- limma_toptable[limma_toptable$coef == '7vs0',]
b<- ifelse(degenes$logFC<0 & degenes$adj.P.Val<0.01, 1, 0)
degenes<- limma_toptable[limma_toptable$coef == '24vs0',]
c<- ifelse(degenes$logFC<0 & degenes$adj.P.Val<0.01, 1, 0)

d<- list()
d$table<- table(a, b, c)
d$labels<- c("2 vs 0","7 vs 0","24 vs 0")


vc<- vennCounts(cbind('2vs0'= a, '7vs0'= b, '24vs0'= c))
vc_perc<- vc[,4] / sum(vc[,4]) * 100
vc<- cbind(as.matrix(vc), vc_perc)

#:> vc
#     2vs0 7vs0 24vs0 Counts    vc_perc
#[1,]    0    0     0  18402 76.2840443
#[2,]    0    0     1   1219  5.0532687
#[3,]    0    1     0   1696  7.0306347
#[4,]    0    1     1   1990  8.2493886
#[5,]    1    0     0    250  1.0363553
#[6,]    1    0     1     32  0.1326535
#[7,]    1    1     0    291  1.2063176
#[8,]    1    1     1    243  1.0073374
#:> 

## Plot:
windows(width= 22/2.54, height= 10/2.54)
par(mfrow= c(1,2), mar= c(0, 0, 1, 2))
plot.venn.diagram(D);  title(main= 'Upregulated')
plot.venn.diagram(d); title(main= 'Downregulated')
savePlot('M:/Documents/LabBook/LabBook_Figures/20110412_vennDiagrams.emf', 'emf') ## See also Powerpoint slide

# ------------------------------------------------------------------------------
# Plot frequency and density of Fold Changes
# ------------------------------------------------------------------------------

# ------------------------------[ Frequency ]-----------------------------------
dev.off()
windows(width= 18/2.54, height= 6/2.54)
par(cex= 0.85, mgp= c(2, 0.75, 0), mfrow= c(1,3), mar= c(4,3,2,0), las= 1, oma=c(0,1,2,0.1), xpd= FALSE)
coefs<- c('2vs0', '7vs0', '24vs0')
clabels<- c('2 hours', '7 hours', '24 hours')
for (i in 1:length(coefs)){
    coef<- coefs[i]
    clab<- clabels[i]
    h<- hist(abs(limma_toptable$logFC[limma_toptable$logFC < 0 & limma_toptable$adj.P.Val < 0.01 & limma_toptable$coef == coef]),
             plot= FALSE, breaks= 50)
    plot(
        smooth.spline(x= h$mids, y= h$counts, spar= 0.45),
        col= 'steelblue', lwd= 2, main= '', xlab= 'Fold change',
        ylab= '', xaxt= 'n', yaxt= 's', xlim= c(0, 5), type= 'l'
    )
    xax<- axis(side= 1, at= 0:5, labels= 0:5)
    h2<- hist(abs(limma_toptable$logFC[limma_toptable$logFC > 0 & limma_toptable$adj.P.Val < 0.01 & limma_toptable$coef == coef]),
        plot= FALSE, breaks= 50)
    lines(
       smooth.spline(x= h2$mids, y= h2$counts, spar= 0.45),
       col= 'firebrick4', lwd= 2
    )
    plot_banner(clab, cex.banner= 1.8, cex.title= 0.85, col.banner= 'grey85', line= 0.15)
}
mtext(side= 3, text= 'Distribution of fold change in LPS-induced and LPS-repressed probes', outer= TRUE, line= 0.5, cex= 0.85, font= 2)
par(mfg=c(1,1))
legend('topright', legend= c('Induced', 'Repressed'), col= c('firebrick4', 'steelblue'), lwd= 2, cex= 0.85, bg= 'white')
mtext(side= 2, text= 'Frequency', line= 2.5, cex= 0.85, las= 0, xpd= NA)
savePlot('M:/Documents/LabBook/LabBook_Figures/20110412_foldchange_frequency.emf', 'emf')

# ------------------------------[ Density ]------------------------------------
dev.off()
windows(width= 18/2.54, height= 6/2.54)
par(cex= 0.85, mgp= c(2, 0.75, 0), mfrow= c(1,3), mar= c(4,0.5,2,0), las= 1, oma=c(0,3,2,0.1), xpd= FALSE)
for (coef in c('2vs0', '7vs0', '24vs0')){
    plot(
        density(abs(limma_toptable$logFC[limma_toptable$logFC < 0 & limma_toptable$adj.P.Val < 0.01 & limma_toptable$coef == coef])),
        col= 'steelblue', lwd= 2, main= '', xlab= 'Fold change', cex.main= 0.85,
        ylab= '', xaxt= 'n', yaxt= 'n', ylim= c(0, 1.6), xlim= c(0, 5)
    )
    yax<-axis(side=2, labels= FALSE, tick= FALSE, at= seq(0, 2, by= 0.5)) 
    xax<- axis(side= 1, at= 0:5, labels= 0:5)
    abline(v= xax, col= 'grey80', lty= 'dotted')
    abline(h= yax, col= 'grey80', lty= 'dotted')
    lines(
       density(abs(limma_toptable$logFC[limma_toptable$logFC > 0 & limma_toptable$adj.P.Val < 0.01 & limma_toptable$coef == coef])),
       col= 'firebrick4', lwd= 2
    )
    plot_banner(coef, cex.banner= 1.8, cex.title= 0.85, col.banner= 'grey85')
    ## mtext(side= 3, text= coef, font= 2, line= 0.5, cex= 0.85, adj= 0)
}
mtext(side= 3, text= 'Distribution of fold change in LPS-induced and LPS-repressed probes', outer= TRUE, line= 0.5, cex= 0.85, font= 2)
par(mfg=c(1,1))
legend('topright', legend= c('Induced', 'Repressed'), col= c('firebrick4', 'steelblue'), lwd= 2, cex= 0.85, bg= 'white')
yax<-axis(side=2, labels= seq(0, 1.5, by= 0.5), at= seq(0, 1.5, by= 0.5), xpd= NA) 
mtext(side= 2, text= 'Density', line= 2.5, cex= 0.85, las= 0, xpd= NA)
savePlot('M:/Documents/LabBook/LabBook_Figures/20110412_foldchange_density.emf', 'emf')


# -----------------------------------------------------------------------------
# USEFUL FUNCTIONS
# -----------------------------------------------------------------------------

plot_banner<- function(plot_title, cex.banner= 2.5, cex.title= par('cex'),
    col.banner= 'grey80', col.border= par('fg'), adj= 0.025, line= 0.1, col.text= par('fg')){
    ## Print on the current plot a banner with a title for the plot 
    ## plot_title= pileup.names[i]
    ## ARGS:
    ## plot_title
    ##     Title to plot
    ## cex.banner
    ##     Expansion character for the height of the banner
    ## cex.title
    ##      Expansion chararacter for the text
    ## col.banner
    ##      Background colour for the banner
    ## col.border
    ##      Colour for the border of the banner
    ## adj
    ##      Horizontal adjustment for the text
    ## line
    ##      Distance of the text from the top border of the graph.
    ## col.text
    ##      Colour for the text
    rect(xleft= par('usr')[1], ybottom= par('usr')[4], 
         xright= par('usr')[2], ytop= par('usr')[4] + strheight('Title', cex= cex.banner, font= 2),
         lwd= 1, border= col.border, col= col.banner, xpd= TRUE)
    mtext(side= 3, line= line, adj= adj, text= plot_title, font= 2, cex= cex.title, col= col.text)
}


# -----------------------------------------------------------------------------
# Script for plotting Venn diagram
# see http://tolstoy.newcastle.edu.au/R/help/03a/1115.html

venn.overlap <-
function(r, a, b, target = 0)
{
#
# calculate the overlap area for circles of radius a and b
# with centers separated by r
# target is included for the root finding code
#
        pi = acos(-1)
        if(r >= a + b) {
                return( - target)
        }
        if(r < a - b) {
                return(pi * b * b - target)
        }
        if(r < b - a) {
                return(pi * a * a - target)
        }
        s = (a + b + r)/2
        triangle.area = sqrt(s * (s - a) * (s - b) * (s - r))
        h = (2 * triangle.area)/r
        aa = 2 * atan(sqrt(((s - r) * (s - a))/(s * (s - b))))
        ab = 2 * atan(sqrt(((s - r) * (s - b))/(s * (s - a))))
        sector.area = aa * (a * a) + ab * (b * b)
        overlap = sector.area - 2 * triangle.area
        return(overlap - target)
}

plot.venn.diagram <-
function(d)
{
#
# Draw Venn diagrams with proportional overlaps
# d$table = 3 way table of overlaps
# d$labels = array of character string to use as labels
#
pi = acos(-1)
csz = 0.1
# Normalize the data
n = length(dim(d$table))
c1 = vector(length = n)
c1[1] = sum(d$table[2, , ])
c1[2] = sum(d$table[, 2, ])
c1[3] = sum(d$table[, , 2])
n1 = c1
#
c2 = matrix(nrow = n, ncol = n, 0)
c2[1, 2] = sum(d$table[2, 2, ])
c2[2, 1] = c2[1, 2]
c2[1, 3] = sum(d$table[2, , 2])
c2[3, 1] = c2[1, 3]
c2[2, 3] = sum(d$table[, 2, 2])
c2[3, 2] = c2[2, 3]
n2 = c2
#
c3 = d$table[2, 2, 2]
n3 = c3
c2 = c2/sum(c1)
c3 = c3/sum(c1)
c1 = c1/sum(c1)
n = length(c1)
# Radii are set so the area is proporitional to number of counts
pi = acos(-1)
r = sqrt(c1/pi)
r12 = uniroot(venn.overlap, interval = c(max(r[1] - r[2], r[2] - r[1],
    0) + 0.01, r[1] + r[2] - 0.01), a = r[1], b = r[2], target = c2[1, 2])$root
r13 = uniroot(venn.overlap, interval = c(max(r[1] - r[3], r[3] - r[1], 0) + 0.01,
    r[1] + r[3] - 0.01), a = r[1], b = r[3], target = c2[1, 3])$root
r23 = uniroot(venn.overlap, interval = c(max(r[2] - r[3], r[3] - r[2], 0) + 0.01,
    r[2] + r[3] - 0.01), a = r[2], b = r[3], target = c2[2, 3])$root
s = (r12 + r13 + r23)/2
x = vector()
y = vector()
x[1] = 0
y[1] = 0
x[2] = r12
y[2] = 0
angle = 2 * atan(sqrt(((s - r12) * (s - r13))/(s * (s - r13))))
x[3] = r13 * cos(angle)
y[3] = r13 * sin(angle)
xc = cos(seq(from = 0, to = 2 * pi, by = 0.01))
yc = sin(seq(from = 0, to = 2 * pi, by = 0.01))
cmx = sum(x * c1)
cmy = sum(y * c1)
x = x - cmx
y = y - cmy
rp=sqrt(x*x + y*y)
frame()
par(usr = c(-1, 1, -1, 1), pty = "s")
## box()
for(i in 1:3) {
    linecol= c('firebrick4', 'dodgerblue', 'seagreen4')
    lines(xc * r[i] + x[i], yc * r[i] + y[i], col= linecol[i], lwd= 2)
}
xl = (rp[1] + (0.7 * r[1])) * x[1]/rp[1]
yl = (rp[1] + (0.7 * r[1])) * y[1]/rp[1]
# text(xl, yl, d$labels[1])                  ## Text group name
# text(xl, yl - csz, d$table[2, 1, 1])     ## Text number 
xl = (rp[2] + (0.7 * r[2])) * x[2]/rp[2]
yl = (rp[2] + (0.7 * r[2])) * y[2]/rp[2]   ## Text number 
# text(xl, yl, d$labels[2])                  ## Text group name
# text(xl, yl - csz, d$table[1, 2, 1])
xl = (rp[3] + (0.7 * r[3])) * x[3]/rp[3]
yl = (rp[3] + (0.7 * r[3])) * y[3]/rp[3]
# text(xl, yl, d$labels[3])                  ## Text group name
# text(xl, yl - csz, d$table[1, 1, 2])
#
# text((x[1] + x[2])/2, (y[1] + y[2])/2, d$table[2, 2, 1]) ## Text number
# text((x[1] + x[3])/2, (y[1] + y[3])/2, d$table[2, 1, 2]) ## Text number
# text((x[2] + x[3])/2, (y[2] + y[3])/2, d$table[1, 2, 2]) ## Text number
#
# text(0, 0, n3)  ## Text common to all
list(r = r, x = x, y = y, dist = c(r12, r13, r23), count1 = c1, count2 =
c2, labels = d$labels)

## Group names
legend( "topright", legend= d$labels, col= linecol, bty= 'n', pch= 1,
    pt.cex= 1.5, pt.lwd= 2 )
}


# -----------------------------------------------------------------------------
# TRITUME
# -----------------------------------------------------------------------------


pig_id <- factor(targets_ori$pig)
time_point <- factor(paste('tp', targets_ori$time_point, sep='.'))

design <- model.matrix(~0 + time_point + pig_id)
fit <- lmFit(limma_affy, design)
fit <- eBayes(fit)
dim(topTable(fit, coef= 6, p.value= 0.01, number= 100000))

> SibShip <- factor(targets$SibShip)
> Treat <- factor(targets$Treatment, levels=c("C","T"))
> design <- model.matrix(~SibShip+Treat)
> fit <- lmFit(eset, design)
> fit <- eBayes(fit)
> topTable(fit, coef="TreatT")


## Select targets at time points:
for (i in c(2, 7, 24)){
    time_1<- 0
    time_2<- i
    
    tp_compare<- which(targets_ori$time_point == time_1 | targets_ori$time_point == time_2)
    targets<- targets_ori[tp_compare, ]
    
    # Subset to have only required timne points
    limma_affy_sub<- limma_affy[, tp_compare] 
    colnames(limma_affy_sub)<- paste(targets$pig, targets$time_point, sep='_')
    limma_affy_sub[1:10,]
    
    # Specify design matrix (see page 38 'Limma user's guide', para '8.3 Paired Samples')
    pig<- factor(targets$pig)
    time_p<- factor(targets$time_point, levels= unique(targets$time_point))
    
    design<- model.matrix(~ 1 + pig + time_p)
    fit <- lmFit(limma_affy_sub, design)
    fit <- eBayes(fit)
    
    time_p_ind<- grep( 'time_p', names(data.frame(fit$coefficients[1:10,])) ) ## Index of the time point tested
    time_p_name<- names(data.frame(fit$coefficients[1:10,]))[time_p_ind]
    
    topTable(fit, coef= time_p_name)  ##
    de_table<- topTable(fit, coef= time_p_name, number= nrow(fit$coefficients), adjust="BH")
    
    comparison<- paste(time_1, 'vs', time_2) ## Make a nice name for this comparison
    de_table$comparison<- comparison
    
    # Write out to file: Data and graphs
    write.table(de_table, 'lwl_bmdm_lps_limma_toptags.txt', sep='\t', row.names=F, append=T, col.names=F)
    
    tiff(paste(time_p_name, '.tiff', sep=''), res= 150, pointsize= 12, units= 'cm', width = 16, height = 9 )
    par(mfrow=c(1,2), cex=0.8, mar= c(5,4,0,1), oma= c(0,0,4,0))
      hist(de_table$P.Value, main= '', xlab='p-value')
      plot(de_table$AveExpr, de_table$logFC, , col= ifelse(de_table$adj.P.Val<0.01, 'blue', 'black'), pch=16, cex=0.2, 
        main= '', xlab= 'Log average expression', ylab= 'log2 fold change')
      title(outer= T, main= paste("Expression differences at", comparison) )
    dev.off()
}


lev <- c(0,2,7,24)
f<- factor(targets_ori$time_point, levels=lev)
design <- model.matrix(~0+f)
colnames(design) <- paste('tp', lev, sep= '.')
fit <- lmFit(limma_affy, design)

cont.tp <- makeContrasts(
    tp.0-tp.2,
    tp.0-tp.7,
    tp.0-tp.24,
    levels=design)
fit2 <- contrasts.fit(fit, cont.tp)
fit2 <- eBayes(fit2)


dim(topTable(fit2, adjust="BH", coef= 1, number= 100000, p.value= 0.01)) ## 854
dim(topTable(fit2, adjust="BH", coef= 2, number= 100000, p.value= 0.01)) ## 4762
dim(topTable(fit2, adjust="BH", coef= 3, number= 100000, p.value= 0.01)) ## 3971

###