# ------------------------------------------------------------------------------
# Comparison between RNSeq and Affymetrix arrays - BMDM samples
# ------------------------------------------------------------------------------

## The BMDM samples to RNAseq is LWL1 see labbook 14/06/2010, time points 0
## and 7h

library(RODBC)
library(biomaRt)
conn<- odbcConnect(dsn= 'pgVitelleschi')


#------------------[ Get probes with an associated transcript ID ]--------------

mart<- useDataset(dataset= "sscrofa_gene_ensembl", useMart("ensembl"))

## All affy probes
affyids<- getBM(attributes= 'affy_porcine', mart= mart)    ## 7443 probes retrieved
dim(affyids) 

probe2transcript<- getBM(
             attributes= c('affy_porcine', 
                           'ensembl_transcript_id'),
             filters= 'affy_porcine',
             value= affyids$affy_porcine,
             mart= mart)
dim(probe2transcript)     ## 8702 probes -> transcript
head(probe2transcript)

## Upload to postgres and generate a table with
## transcript_id, ctrl/lps, probes exprs, FPKM
sqlSave(conn, probe2transcript, rownames= FALSE)

## Table matching RNAseq with Affymetrix
gene_expr<- sqlQuery(conn, "
    select transcript_id, time_point, log2_intensity, fpkm
    from cufflinks_transcript_gtf inner join 
        -- Extract probes with a transcript_id get only 0 and 7 hours and only pig LWL1
        ( 
        select * from probe2transcript inner join affyarray_bmdm_avg on
             probe_set_id = affy_porcine 
        where time_point in (0, 7) and pig = 'LWL1'
        ) as affy on 
        -- Link probe to transcripts
        transcript_id = ensembl_transcript_id AND
        -- Link time points in affy and cufflinks
        time_point = (CASE WHEN source = '20110207_bmdm_ctrl' THEN 0 WHEN source = '20110207_bmdm_lps' THEN 7 END)::int
    where 
        -- Filter cufflinks table
        source in ('20110207_bmdm_ctrl', '20110207_bmdm_lps') AND feature like 'transcript'
")
head(gene_expr)



# ------------------------------------------------------------------------------
# Plot comparison
# ------------------------------------------------------------------------------

rnaseq_ctrl<- gene_expr[gene_expr$time_point == 0, 'fpkm']
affy_ctrl<- gene_expr[gene_expr$time_point == 0, 'log2_intensity']
rnaseq_lps<- gene_expr[gene_expr$time_point == 7, 'fpkm']
affy_lps<- gene_expr[gene_expr$time_point == 7, 'log2_intensity']
head(gene_expr)

library("geneplotter")  ## from BioConductor
require("RColorBrewer") ## from CRAN

source('U:/Documents/ScriptArchive/R/plot_template_single.R')

## LPS
x<- affy_lps[rnaseq_lps > 0]
y<- log2(rnaseq_lps[rnaseq_lps > 0])
mylm<- lm(y ~ x + I(x^2) + I(x^3))
summary(mylm)
newx<- seq(2, 20) 
prd<- predict(mylm, newdata= data.frame(x=newx), interval = c("prediction"), type="response")

tiff('M:/Documents/LabBook/LabBook_Figures/20110331_affy_vs_rnaseq_lps.tif', res= 150, pointsize= 8, units= 'cm',
    width = 7, height = 8)
smoothScatter(affy_lps, log2(rnaseq_lps + 0.01), bandwidth=0.01,
    xlab= 'log2(Affymetrix)', ylab= 'log2(RNAseq)', main= 'Affymetrix vs RNAseq (7 hours)')
lines(newx, prd[,1], col="red", lty=2, lwd= 1.5)
dev.off()


## CTRL
x<- affy_ctrl[rnaseq_ctrl > 0]
y<- log2(rnaseq_ctrl[rnaseq_ctrl > 0])
mylm<- lm(y ~ x + I(x^2) + I(x^3))
summary(mylm)
newx<- seq(2, 20) 
prd<- predict(mylm, newdata= data.frame(x=newx), interval = c("prediction"), type="response")

tiff('M:/Documents/LabBook/LabBook_Figures/20110331_affy_vs_rnaseq_ctrl.tif', res= 150, pointsize= 8, units= 'cm',
    width = 7, height = 8)
smoothScatter(affy_ctrl, log2(rnaseq_ctrl + 0.01), bandwidth=0.01,
    xlab= 'log2(Affymetrix)', ylab= 'log2(RNAseq)', main= 'Affymetrix vs RNAseq (0 hours)')
lines(newx, prd[,1], col="red", lty=2, lwd= 1.5)
dev.off()

    
cor.test( x= log2(rnaseq_ctrl + 0.01 ), y= affy_ctrl )
cor.test( x= log2(rnaseq_lps + 0.01  ), y= affy_lps)

cor( cbind(rnaseq_ctrl, rnaseq_lps, affy_ctrl, affy_lps) )
cor( cbind(rnaseq_ctrlrnaseq_ctrl, rnaseq_lps, affy_ctrl, affy_lps) )

sqlDrop(conn, probe2transcript)

# ------------------------------------------------------------------------------
# Get cDNA sequence of RNAseq expressed genes
# ------------------------------------------------------------------------------

rnaseq_transcripts<- sqlQuery(conn, "
    select distinct transcript_id from cufflinks_transcript_gtf
    where source in ('20110207_bmdm_ctrl', '20110207_bmdm_lps') AND
          feature like 'transcript' AND
          fpkm > 0
", stringsAsFactors= FALSE)[,1]
rnaseq_transcripts[1:10]

writeLines(rnaseq_transcripts[,1], 'U:/Tritume/rnaseq_transcripts.txt')
## rnaseq_transcripts<- readLines('/mnt/ris-fas1a/users/dberaldi/Tritume/rnaseq_transcripts.txt')
mart<- useDataset(dataset= "sscrofa_gene_ensembl", useMart("ensembl"))

for(i in seq(1, length(rnaseq_transcripts), by= 1000)){
## This is quite stupid: Retrieve sequences in chunks of 1000 otherwise
## all at once biomart chokes (also the web version).
    from<- i
    to<- i + 999
    if(to > length(rnaseq_transcripts)){
        to<- length(rnaseq_transcripts)
    }
    cdna<- getBM(	
             attributes= c("cdna", "ensembl_transcript_id"), 
             filters= "ensembl_transcript_id",
             value= rnaseq_transcripts[from:to],
             mart= mart)
    write.table(cdna, 'D:/Tritume/rnaseq_transcript_sequences.txt', col.names= FALSE, row.names= FALSE, quote= FALSE, append= TRUE)
}
cdna[1:2,1:2]

# ------------------------------[ TRITUME ]-------------------------------------

## Retrieve affymetrix expression from postgres
affy<- sqlQuery(conn, "
    SELECT * FROM affyarray_bmdm_avg INNER JOIN 
    WHERE pig = 'LWL1' AND
        time_point in (0, 7)
")
head(affy)
## Get preobes with associated transcript_id
affy_transcript
affy_transcript<- affy[affy$probe_set_id %in% probe2transcript$affy_porcine, ]

head(affy_transcript)
dim(affy_transcript)  ## 14886



# ------------------------------------------------------------------------------
# Get transcript expression by RNAseq (0 and 7 h)
# ------------------------------------------------------------------------------

rnaseq<- sqlQuery(conn, "
    SELECT * FROM cufflinks_transcript_gtf
    WHERE source in ('20110207_bmdm_ctrl', '20110207_bmdm_lps') AND
          feature like 'transcript';
")
head(rnaseq)
