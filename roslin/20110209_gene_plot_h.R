# -----------------------------------------------------------------------------
# Plot CAGE tags around the inferred TSSs
# -----------------------------------------------------------------------------

library(RODBC)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Mmusculus.UCSC.mm9")
library("BSgenome.Sscrofa.UCSC.susScr2")

conn<- odbcConnect(dsn= 'pgVitelleschi')

setwd('M:/Documents/LabBook/LabBook_Figures/20110209_immunogenes')


# -------------------------------[ CSF1R humans ]------------------------------

gene_name<- 'CSF1R'
species<- 'hsapiens'
chr<- '5'
start<- 149434383 - 1000
end<- 149492935 + 1000
strand<- '-'
## Get sequence
genome<- Hsapiens

csf1r_hs<- sqlQuery(conn,
  paste("select * from tss_bed where chr = '", chr, "' and strand = '", strand, "' AND tss_pos between ", start, " and ", end , " AND species = '", species, "'", sep='')
    )
transcript_starts<- sqlQuery(conn,
    paste("select * from transcript_attributes where reference_gene_name = '", gene_name, "' and species = '", species, "'", sep= '')
    )
plot_genes()
# -------------------------------[ CSF1R mouse ]-------------------------------

gene_name<- 'CSF1R'
species<- 'mmusculus'
chr<- 18
start<- 61265226 - 1000
end<- 61265522 + 1000
strand<- '+'
genome<- Mmusculus

csf1r_hs<- sqlQuery(conn,
  paste("select * from tss_bed where chr = '", chr, "' and strand = '", strand, "' AND tss_pos between ", start, " and ", end , " AND species = '", species, "'", sep='')
    )
transcript_starts<- sqlQuery(conn,
    paste("select * from transcript_attributes where reference_gene_name = '", gene_name, "' and species = '", species, "'", sep= '')
    )
plot_genes()
# -------------------------------[ CSF1R pig ]---------------------------------
gene_name<- 'CSF1R'
species<- 'sscrofa'
chr<- 2
start<- 136552808-1000
end<- 136552808+1000
strand<- '-'
genome<- Sscrofa

csf1r_hs<- sqlQuery(conn,
  paste("select * from tss_bed where chr = '", chr, "' and strand = '", strand, "' AND tss_pos between ", start, " and ", end , " AND species = '", species, "'", sep='')
    )
transcript_starts<- sqlQuery(conn,
    paste("select * from transcript_attributes where reference_gene_name = '", gene_name, "' and species = '", species, "'", sep= '')
    )

plot_genes()

# -------------------------------[ SPI1 humans ]------------------------------

gene_name<- 'SPI1'
species<- 'hsapiens'
chr<- '11'
start<- 47400127 - 1000
end<- 47400127 + 1000
strand<- '-'
## Get sequence
genome<- Hsapiens

csf1r_hs<- sqlQuery(conn,
  paste("select * from tss_bed where chr = '", chr, "' and strand = '", strand, "' AND tss_pos between ", start, " and ", end , " AND species = '", species, "'", sep='')
    )
transcript_starts<- sqlQuery(conn,
    paste("select * from transcript_attributes where reference_gene_name = '", gene_name, "' and species = '", species, "'", sep= '')
    )

plot_genes()
# -------------------------------[ SPI1 mouse ]-------------------------------

gene_name<- 'SPI1'
species<- 'mmusculus'
chr<- '2'
start<- 90936253 - 1000
end<- 90936937 + 1000
strand<- '+'
genome<- Mmusculus

csf1r_hs<- sqlQuery(conn,
  paste("select * from tss_bed where chr = '", chr, "' and strand = '", strand, "' AND tss_pos between ", start, " and ", end , " AND species = '", species, "'", sep='')
    )
transcript_starts<- sqlQuery(conn,
    paste("select * from transcript_attributes where reference_gene_name = '", gene_name, "' and species = '", species, "'", sep= '')
    )

plot_genes()
# -------------------------------[ SPI1 pig ]---------------------------------

gene_name<- 'SPI1'
species<- 'sscrofa'
chr<- 2
start<- 13254376 - 1000
end<- 13254376 + 1000
strand<- '-'
genome<- Sscrofa

csf1r_hs<- sqlQuery(conn,
  paste("select * from tss_bed where chr = '", chr, "' and strand = '", strand, "' AND tss_pos between ", start, " and ", end , " AND species = '", species, "'", sep='')
    )
transcript_starts<- sqlQuery(conn,
    paste("select * from transcript_attributes where reference_gene_name = '", gene_name, "' and species = '", species, "'", sep= '')
    )

plot_genes()


# -------------------------------[ TNF humans ]------------------------------

gene_name<- 'TNF'
species<- 'hsapiens'
chr<- '6'
start<- 31543003 - 1000
end<- 31543344 + 1000
strand<- '+'
## Get sequence
genome<- Hsapiens

csf1r_hs<- sqlQuery(conn,
  paste("select * from tss_bed where chr = '", chr, "' and strand = '", strand, "' AND tss_pos between ", start, " and ", end , " AND species = '", species, "'", sep='')
    )
transcript_starts<- sqlQuery(conn,
    paste("select * from transcript_attributes where reference_gene_name = '", gene_name, "' and species = '", species, "'", sep= '')
    )

plot_genes()
# -------------------------------[ TNF mouse ]-------------------------------

gene_name<- 'TNF'
species<- 'mmusculus'
chr<- '17'
start<- 35338952 - 1000
end<- 35338952 + 1000
strand<- '-'
genome<- Mmusculus

csf1r_hs<- sqlQuery(conn,
  paste("select * from tss_bed where chr = '", chr, "' and strand = '", strand, "' AND tss_pos between ", start, " and ", end , " AND species = '", species, "'", sep='')
    )
transcript_starts<- sqlQuery(conn,
    paste("select * from transcript_attributes where reference_gene_name = '", gene_name, "' and species = '", species, "'", sep= '')
    )

plot_genes()
# -------------------------------[ TNF pig ]---------------------------------

gene_name<- 'TNF'
species<- 'sscrofa'
chr<- 7
start<- 27658519 - 1000
end<- 27658519 + 1000
strand<- '+'
genome<- Sscrofa

csf1r_hs<- sqlQuery(conn,
  paste("select * from tss_bed where chr = '", chr, "' and strand = '", strand, "' AND tss_pos between ", start, " and ", end , " AND species = '", species, "'", sep='')
    )
transcript_starts<- sqlQuery(conn,
    paste("select * from transcript_attributes where reference_gene_name = '", gene_name, "' and species = '", species, "'", sep= '')
    )

plot_genes()



# -----------------------------------------------------------------------------
# Plot broad region
# -----------------------------------------------------------------------------
plot_genes<- function(){
    windows(width= 24/2.54, height= 6/2.54)
    par(mar= c(0,0,0,1), mfrow= c(1,2), oma=c(2,2,0,0), bg= 'white' )
    
    plot(x=csf1r_hs$tss_pos, y=csf1r_hs$tpm, type='h', xlab= '', ylab='', main= '', xlim=c(start, end),
         frame.plot= FALSE, cex.axis= 0.65, xaxt= 'n', lwd= 2, yaxt= 'n')
    #at<- seq(start, end, length.out= 11)
    x<- axis(side= 1, labels= NA)
    mtext(side= 1, line= 0.25, at= x, text= x, cex= 0.65)
    y<- axis(side= 2, labels= NA)
    mtext(side= 2, line= 0.55, at= y, text= y, cex= 0.65, las= 1)
    
#    mtext(side= 1, line= 2.5, text= 'Chromosome position', cex= 0.65)
#    mtext(side= 2, line= 2.5, text= 'Tags per million', cex= 0.65)
    abline(v= transcript_starts$transcript_5p_end, col= 'blue', lty= 'dotted')
 
    legend('topleft', legend= paste(gene_name, ' - ', species, sep= ''), text.col= 'firebrick4', bg= 'white', box.col= 'white', inset= c(0.01, 0))

     
    # ---------------------------[ Narrow around peak ]------------------------

    peak_tpm<- max(csf1r_hs$tpm) ## Highest peak
    peak_pos<- csf1r_hs$tss_pos[ which(csf1r_hs$tpm == peak_tpm) ] ## Position of peak
    
    ## Offset around the peak
    peak_left<- peak_pos - 40
    peak_right<- peak_pos + 40
    
    ## Subset to have only narrow region. Add 0 to positions with no tags
    csf1r_peak_hs<- csf1r_hs[csf1r_hs$tss_pos > peak_left & csf1r_hs$tss_pos < peak_right, ]
    pos<- seq(peak_left, peak_right)
    tpm<- rep(0, length(pos))
    tpm[pos %in% csf1r_peak_hs$tss_pos]<- csf1r_peak_hs$tpm
    
    tss<- strsplit(getSeq(genome, names= paste('chr', chr, sep=''), start= peak_left, end= peak_right, strand= strand), split= '')[[1]]
    if (strand == '-') ## Reverse the tpm vetcor if it comes from (-)strand
        tpm<- rev(tpm)
    
    col<- tss
    col[tss == 'A'] <- 'green'
    col[tss == 'C'] <- 'orange'
    col[tss == 'G'] <- 'red'
    col[tss == 'T'] <- 'blue'
    
    plot(x=pos, y= tpm, type='h', xlab= '', lwd= 2,
         xaxt= 'n', yaxt= 'n', ylab='Tags per million', main= '', frame.plot= FALSE, col= col, cex.axis= 0.65, cex.lab= 0.65)
    mtext(text= tss, at= seq(peak_left, peak_right), side= 1, cex= 0.6, line= -0.5, col= col)
    y<- axis(side= 2, labels= NA)
    mtext(side= 2, line= 0.55, at= y, text= y, cex= 0.65, las= 1)
 
    savePlot(paste(gene_name, '_', species, '.emf', sep= ''), 'emf')
    dev.off()
    }

# ------------------------------[ TRITUME ]------------------------------------
