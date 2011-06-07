#
#  Sequences of TSSs for various genes in pig/mouse/human
#  Output in a FASTA file
#  These sequences used for multialign

library(RODBC)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Mmusculus.UCSC.mm9")
library("BSgenome.Sscrofa.UCSC.susScr2")

setwd('M:/Documents/LabBook/LabBook_Figures/20110209_immunogenes')

conn<- odbcConnect(dsn= 'pgVitelleschi')

genomes<- c(Hsapiens, Mmusculus, Sscrofa)
species<- c('hsapiens', 'mmusculus', 'sscrofa')
gene_names<- c('CSF1R', 'SPI1', 'TNF')

tss_genes<- sqlQuery(conn, 'select * from tss_genes');

## Where fasta file will be
fasta_file<- 'tss_sequences.fa'
file.remove(fasta_file)

for (gene_name in gene_names){
    gene<- tss_genes[tss_genes$reference_gene_name == gene_name, ]
    for ( i in 1:length(species) ){
        chr<- gene$chr[gene$species == species[i] ]
        peak<- gene$tss_pos[gene$species == species[i] ]
        peak_left<- peak - 40
        peak_right<- peak + 40
        strand<- as.character(gene$strand[gene$species == species[i] ])
        sequence<- getSeq(genomes[i][[1]], names= paste('chr', chr, sep=''), start= peak_left, end= peak_right, strand= strand)
        write.table( file= fasta_file, row.names= FALSE, col.names=FALSE, quote=FALSE, append= TRUE,
            x= paste('>', gene_name, '_', species[i], '_', chr, ':', peak_left, '-', peak_right, ':', strand, '_peak_at:', peak ,sep='')
            )
        write.table( file= 'D:/Tritume/tss_sequences.fa', row.names= FALSE, col.names=FALSE, quote=FALSE, append= TRUE,
            x= sequence
            )
        }
    ## Add a blank line between genes
    write.table( file= 'D:/Tritume/tss_sequences.fa', row.names= FALSE, col.names=FALSE, quote=FALSE, append= TRUE, x= '' )
}