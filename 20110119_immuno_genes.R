library(RODBC)
conn<- odbcConnect(dsn= 'pgVitelleschi')

# ------------------------------------------------------------------------------
# Plot TSS peaks for each gene of interest.
# ------------------------------------------------------------------------------

tss_ig<- sqlFetch(conn, 'tss_peaks_immunogenes')
names(tss_ig)
gene<- 'IDO1' ## ACP5  ACTB  CAPG  CCL20 CEBPB CSF1R GSN   IDO1  IGF1  IL12A IL12B IL1B  IL6   MITF  SPI1  STAT4 TNF

genes<- unique(tss_ig$reference_gene_name)
xlim= c(-200,200)
n<- 1
windows(width= 23/2.54, height= 14/2.54, bg= 'white')
par(mfcol= c(3,3), mar= c(2,2,1,1), oma= c(4,2,2,0))
for ( gene in genes ){
    for (species in c('hsapiens', 'mmusculus', 'sscrofa' ) ){
        tpm<- tss_ig[tss_ig$reference_gene_name == gene & tss_ig$species == species, 'tpm']
        pos<- tss_ig[tss_ig$reference_gene_name == gene & tss_ig$species == species, 'tss_reftrans_dist']
        if( is.na(tpm[1]) ){ tpm<- rep(0, ((1+xlim[2]) - xlim[1]) ); pos<- xlim[1]:xlim[2] }

        plot(type= 'n', x= pos, y= tpm, ylim= c(0, max(tpm)+1), ylab='', xlim= xlim, frame.plot= F)
        abline(v= 0, lty= 'solid', col= 'lightblue', lwd= 4)
        if ( all(tpm == 0) == FALSE )
            points(type= 'h', x= pos, y= tpm, lwd= 2)
        }
    legend('topright', legend=gene, bty= 'n', cex= 1.25)
    mtext(outer= TRUE, "Tags per million (TPM)", side= 2, cex= 0.8, line= 0.2)
    mtext(outer= TRUE, "Position relative to 5' end", side= 1, cex= 0.8, line= +1)
    if (n %% 3 == 0){
        savePlot(paste('M:/Documents/LabBook/LabBook_Figures/20110126_immunogenes/tss_', n, '.emf'), 'emf')        
        }  
    n<- n + 1
    }
savePlot(paste('M:/Documents/LabBook/LabBook_Figures/20110126_immunogenes/tss_', n, '.emf'), 'emf')

# -----------------------------------------------
#  3D plots
# -----------------------------------------------
require(scatterplot3d)

gene<- 'IGF1'

xlim= c(-200, 200)
genes<- unique(tss_ig$reference_gene_name)

windows(width= 23/2.54, height= 20/2.54)
lyout<- layout(matrix(c(1,3,2,4), nrow=2))
layout.show(lyout)
par(xpd= T)
n<- 1
for (gene in genes){
    gene3d<- tss_ig[tss_ig$reference_gene_name == gene,
                     c('tss_reftrans_dist', 'species', 'tpm')]
    ## Add a row with 0 to plot missing data
    for (species in unique(tss_ig$species) ){
        gene3d<- rbind(gene3d, c(0, species, 0))
        }    
    species_col<- as.character(gene3d$species)
    species_col[species_col == 'hsapiens']<- 'dodgerblue'
    species_col[species_col == 'mmusculus']<- 'firebrick'
    species_col[species_col == 'sscrofa']<- 'darkolivegreen4'

    scatterplot3d(gene3d, type="h", lwd= 2, pch=" ", color=species_col, box= F, 
        angle= 65, lty.grid= 'dotted', col.grid= 'lightgrey', scale.y= 0.35, xlim= xlim, zlim= c(0, max(as.numeric(gene3d$tpm)) + 1),
        ylab= '', xlab= '', zlab= '', y.ticklabs= '', cex.lab= 0.95, mar= c(2,2,0,0))
    legend('topleft', legend= gene, bty= 'n', inset= c(0.05, 0.15), cex= 1.5)
#    if (n == 3) {
        legend('topright', legend= c('Pig', 'Mouse', 'Human'), col= c('darkolivegreen4', 'firebrick', 'dodgerblue'), 
            lwd= 2, bty= 'n', cex= 0.95, inset= c(0.05, 0.15)) 
        mtext(text= "Position relative to 5' end", side= 1, line= 0.75, cex= 0.9)
        mtext(text= "Tags per million", side= 2, line= 0.85, cex= 1.1)
 #       }
    if (n == 4) {
        savePlot(paste('M:/Documents/LabBook/LabBook_Figures/20110126_immunogenes/tss_', gene, '_3D.emf'), 'emf')
        n<- 0
        }
    n<- n + 1
    }
savePlot(paste('M:/Documents/LabBook/LabBook_Figures/20110126_immunogenes/tss_', gene, '_3D.emf'), 'emf')

# ------------------------------------------------------------------------------
# Fetch sequences surrounding each peak
# ------------------------------------------------------------------------------    
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Mmusculus.UCSC.mm9")
library("BSgenome.Sscrofa.UCSC.susScr2")

peaks<- sqlFetch(conn, 'tss_genes')
head(peaks)
offsets<- c(-200, 200) ## Bases to fetch from left and right of the peak
tss_seq<- vector() ## Where sequences will be, in order!
gene_id<- vector()
gene_name<- vector()
species_vect<- vector()
tss_id<- vector()
for (species in (unique(peaks$species) ) ){
    if (species == 'hsapiens') genome<- Hsapiens
    if (species == 'mmusculus') genome<- Mmusculus
    if (species == 'sscrofa') genome<- Sscrofa
    message(paste('Species: ', species))
    peak_pos<- peaks[peaks$species == species,]
    for ( i in 1:nrow(peak_pos) ){
        message(paste('Fetching TSS sequence for', peak_pos$ensembl_gene_id[i]))
        peak_seq<- getSeq(genome, names= paste('chr', peak_pos$chr[i], sep= ''), start= peak_pos$tss_pos[i] + offsets[1], end= peak_pos$tss_pos[i] + offsets[2], strand= peak_pos$strand[i])
        tss_seq<- c(tss_seq, peak_seq)
        gene_id<- c(gene_id, as.character(peak_pos$ensembl_gene_id[i]))
        gene_name<- c(gene_name, as.character(peak_pos$reference_gene_name[i]))
        species_vect<- c(species_vect, as.character(peak_pos$species[i]))
        tss_id<- c(tss_id, as.character(peak_pos$tss_id[i]))
        }
    }
tss_sequences<- data.frame(cbind(species= species_vect, ensembl_gene_id= gene_id, reference_gene_name= gene_name, tss_peak_id= tss_id, sequence= tss_seq))
write.table(tss_sequences, 'D:/Tritume/peak_sequences.txt', sep= '\t', row.names= F)
sqlQuery(conn, "select read_table($$ file: 'D:/Tritume/peak_sequences.txt', table: 'peak_sequences', header:True, overwrite:True $$)")
sqlQuery(conn, "comment on table peak_sequences is $$Sequence around the TSS peaks as in tss_genes. +/-200 bp. Sequence is reverse complemented when it comes from the minus strand. See 20110119_immuno_genes.R $$ ")


# -----------------------------
# TRITUME
# -----------------------------

      
tpm<- tss_ig[tss_ig$reference_gene_name == gene, 'tpm']
pos<- tss_ig[tss_ig$reference_gene_name == gene, 'tss_reftrans_dist']
species<- tss_ig[tss_ig$reference_gene_name == gene, 'species']



my.mat <- matrix(runif(25), nrow=5)
  dimnames(my.mat) <- list(LETTERS[1:5], letters[11:15])

  s3d.dat <- data.frame(cols=as.vector(col(my.mat)),
      rows=as.vector(row(my.mat)),
      value=as.vector(my.mat))
  scatterplot3d(s3d.dat, type="h", lwd=5, pch=" ",
      x.ticklabs=colnames(my.mat), y.ticklabs=rownames(my.mat),
      color=grey(25:1/40), main="scatterplot3d - 4")



df1<- data.frame(xv= NA, yv= NA)
df1<- df1[0,]
df1<- rbind(df1, c(1,2), c(3,4))