#!/usr/bin/R

# ----------------------------------[ README ]---------------------------------
# 26/05/2011
# dario.beraldi@roslin.ed.ac.uk
# 
# Analysis of peak mass in chip-seq libraries to detect differential peaks
#
# If DESeq is not installed do:
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("DESeq")
#
# USAGE:
# Change INPUT section below as appropriate. Then do:
# R --slave --no-save --no-restore --no-environ < 20110526_deseq_chipseq_markb.R
# 
# NOTE: Ths script is for two sample only. Some changes (not many) to analyze more
# samples, but two conditions only.
#
# Input files should look like this:
# ------------------------------------------------------------------------------
# chr1	chip	peak	82835671	82838071	.	+	.	 peak_id "Merged-chr1-82836871-2"; file "con_k4b.peaks|lps_k4b.peaks";	1030	2369	2401	0.9866722
# chr1	chip	peak	99614386	99615918	.	+	.	 peak_id "Merged-chr1-99615152-2"; file "con_k4b.peaks|lps_k4b.peaks";	244	1016	1533	0.6627527
# chr1	chip	peak	120586129	120587129	.	+	.	 peak_id "2-chr1-5856"; file "lps_k4b.peaks";	106	958	1001	0.9570429
# chr1	chip	peak	4846590	4850158	.	+	.	 peak_id "Merged-chr1-4848374-2"; file "con_k4b.peaks|lps_k4b.peaks";	3049	3552	3569	0.9952368
# chr1	chip	peak	7077740	7081382	.	+	.	 peak_id "Merged-chr1-7079561-2"; file "con_k4b.peaks|lps_k4b.peaks";	2296	3483	3643	0.9560801
# chr1	chip	peak	10221664	10224222	.	+	.	 peak_id "Merged-chr1-10222943-2"; file "con_k4b.peaks|lps_k4b.peaks";	1964	2519	2559	0.9843689
# ------------------------------------------------------------------------------


# -----------------------------[ INPUT ]----------------------------------------

## These are GTF files typically produced by the pipeline findPeaks > mergePeaks > coverageBed.
con_peaks<- '/exports/work/vet_roslin_nextgen/markb/homer/output/20110519_h3k4me3_mb/con_k4b.coverageBed.gtf'
lps_peaks<- '/exports/work/vet_roslin_nextgen/markb/homer/output/20110519_h3k4me3_mb/lps_k4b.coverageBed.gtf'

## Conditions: A vector of conditions for each sample:
conds<- c('ctrl', 'lps')

## A name for the output file:
outfile<- 'deseq_nbinomtest.txt'

# -------------------[ DESeq Differential precipitation ]----------------------

library(DESeq)

## Read data
cat('Reading datafiles...\n')
con_df<- read.table(con_peaks, header= FALSE, sep= '\t', stringsAsFactor= FALSE, quote= '')
lps_df<- read.table(lps_peaks, header= FALSE, sep= '\t', stringsAsFactor= FALSE, quote= '')

##
attr_start_marker<- 'peak_id "'
attr_end_marker<- '"; file "'
get_peak_id<- function(attribute, attr_start_marker= attr_start_marker, attr_end_marker= attr_end_marker){
    ## Extract the peak_id from gtf attribute line:
    ## attribute<- 'peak_id "Merged-chr1-82836871-2"; file "con_k4b.peaks|lps_k4b.peaks";'
    ## Peak id is flanked by the markers attr_start_marker and attr_end_marker
    plen<- nchar(attr_start_marker)
    id_start<- regexpr(attr_start_marker, attribute)
    id_end<- regexpr(attr_end_marker, attribute)-1
    peak_id<- paste(strsplit(attribute, split= '')[[1]][(id_start+plen):id_end], collapse= '')
    return(peak_id)
}

## Extract peak identifiers
cat('Extracting peak identifiers...\n')
con_df_peaks<- unlist(lapply(con_df$V9, get_peak_id, attr_start_marker, attr_end_marker))
lps_df_peaks<- unlist(lapply(lps_df$V9, get_peak_id, attr_start_marker, attr_end_marker))

if(all(con_df_peaks == lps_df_peaks) == FALSE){
    stop('Peak IDs in the two dataframes are not identical')
} 


# ------------------------------[ Prepare for DESeq ]---------------------------
## Same as edger
cat('Preparing count table...\n')
countsTable<- as.data.frame(cbind(
    con_df$V10,
    lps_df$V10
))
names(countsTable)<- conds
rownames(countsTable)<- con_df_peaks

# Counts table should look like this:
#
cat('Count table:\n\n')
head(countsTable)
#                       ctrl  lps
# Merged-chr1-82836871-2  792 1030
# Merged-chr1-99615152-2  349  244
# 2-chr1-5856              20  106
# Merged-chr1-4848374-2  2414 3049
# Merged-chr1-7079561-2  2325 2296
# Merged-chr1-10222943-2 1890 1964
# Merged-chr1-16646780-2 2321 2262
# Merged-chr1-36569453-2  997  985
# Merged-chr1-37486957-2 3002 3364
# Merged-chr1-52821962-2  612  496

# ------------------------------------------------------------------------------
# DESeq Differential Expression
# ------------------------------------------------------------------------------

cat('\nTesting for differential counts...\n')
cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )       ## Factor to adjust library sizes to account for differences in seq depth

## For each gene, use the mean across all samples to estimate the associated variance
## method= 'blind' allows for no replication.
## pool= TRUE is the same of blind. (pool is deprected in newer versions)
cds <- estimateVarianceFunctions( cds, pool= TRUE )

res <- nbinomTest( cds, unique(conds)[1], unique(conds)[2])  ## Binomial test for D.E.

cat('Differential test table:\n\n')
head(res)

## Send to output
cat('\nWriting out data...\n')
write.table(res, file= outfile, col.names= TRUE, row.names= FALSE, sep= '\t', quote= FALSE)

# :> head(res)
#                      id  baseMean baseMeanA  baseMeanB foldChange log2FoldChange      pval padj resVarA resVarB
#1 chr1_3395081_3396335_+  102.7475  139.1761   66.31884  0.4765101    -1.06942140 0.1710383    1      NA      NA
#2 chr1_4774332_4777028_+ 2072.5642 2180.7832 1964.34522  0.9007522    -0.15079784 0.5572685    1      NA      NA
#3 chr1_4797077_4799629_+ 1622.4547 1794.3017 1450.60776  0.8084525    -0.30676503 0.2705804    1      NA      NA
#4 chr1_4846055_4850645_+ 2738.4182 2596.1704 2880.66603  1.1095828     0.15001731 0.5486946    1      NA      NA
#5 chr1_5008334_5011146_+ 1008.2295  906.7861 1109.67291  1.2237427     0.29130025 0.4413344    1      NA      NA
#6 chr1_5072857_5075547_+ 1359.9277 1339.3028 1380.55266  1.0307995     0.04376377 0.9058069    1      NA      NA
#:> 

# plot( res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( res$padj < .1, "red", "black" ),
#      ylab= 'Log2 fold change', xlab= 'Base mean', cex.main= 0.85, main= 'Peaks found D.E. (padj <.1) in\nrelation to fold change and base mean')
# abline(h= c(-2, 2), col= 'dodgerblue', lwd= 2)
# savePlot('deseq_basemean_diff.bmp', 'bmp')

cat('Done!\n')