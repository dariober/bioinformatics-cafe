#
# Summary stats and graph for SNPs from AM and BMDM
# See also:
#    /exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_xxx for alignments
#    /exports/work/vet_roslin_nextgen/dario/samtools/output/20110414_rnaseq_mpileup
#

library(RODBC)
conn<- odbcConnect(dsn='pgVitelleschi')
odbcClose(conn)

## Number of bases covered by depth of sequences:
## See Labbook 14/04/2011
nreads_am<- c(321029912, 282392270, 217568701, 158747809, 41671388, 19655707)
nreads_bmdm<- c(360720601, 312831535, 229306738, 163027660, 46918762, 20659379)
depth<- rep(c(1, 2, 5, 10, 50, 100))
windows(height= 9/2.54, width= 9/2.54)
par(las= 1, fg= 'grey25', cex= 0.75, mgp=c(1.75, 0.75, 0), mar= c(4,4,2,1))
plot(x= depth, y= nreads_bmdm/10^6, col= 'dodgerblue', pch= 19, type= 'o', ylab= '', main= 'Number of positions by depth')
points(x= depth, y= nreads_am/10^6, col= 'firebrick4', pch= 19, type= 'o')
mtext(side= 2, 'No. positions (millions)', las= 0, cex= 0.75, line= 2.75)
legend('topright', legend= c('Large white', 'Large w. x Landr.'), col= c('dodgerblue', 'firebrick4'), lwd= 1.2 )
savePlot('M:/Documents/LabBook/LabBook_Figures/20110415_rnaseq_positions.emf', 'emf')

## Number of SNPs
snp.all<- sqlQuery(conn, 'SELECT file, qual FROM vcf')
snp<- with(snp.all[snp.all$qual > 0, ], { aggregate(qual, by= list(file), length ) })
snp.10<- with(snp.all[snp.all$qual > 10, ], { aggregate(qual, by= list(file), length ) })
snp.20<- with(snp.all[snp.all$qual > 20, ], { aggregate(qual, by= list(file), length ) })

windows(height= 10/2.54, width= 9/2.54)
par(cex= 0.75, xpd= TRUE)
bp<- barplot(snp$x/1000, ylab= '', names.arg= c('AM', 'ctrl', 'lps', 'BMDM', 'ctrl', 'lps'),
    border= rep(c('steelblue2', 'firebrick2'), each= 3), las= 1,  col=rep(c('steelblue2', 'firebrick2'), each= 3))
mtext(side= 2, text= 'No. SNPs (x1000)', line= 2.75, cex= 0.75)
mtext(side= 3, text=, 'Number of SNPs (phred >3, >10, >20)', font= 2, cex= 1, line= 2.5)
legend('topright', bty= 'n', inset= c(-0.14, -0.1), legend= c('Large white', 'Large w. x Landr.'), col= c('steelblue4', 'firebrick4'), pch= 19 )

barplot(snp.10$x/1000, border= rep(c('steelblue3', 'firebrick3'), each= 3), col=rep(c('steelblue3', 'firebrick3'), each= 3), add= TRUE, yaxt= 'n')
barplot(snp.20$x/1000, border= rep(c('steelblue4', 'firebrick4'), each= 3), col=rep(c('steelblue4', 'firebrick4'), each= 3), add= TRUE, yaxt= 'n')
savePlot('M:/Documents/LabBook/LabBook_Figures/20110415_rnaseq_no_snp.emf', 'emf')
graphics.off()

# -----------------------------------------------------------------------------
# Annotation from ANNOVAR:
# -----------------------------------------------------------------------------
## See Labbook 16/04/2011

snp_ann<- sqlQuery(conn, "
    SELECT file, region, count(distinct chrom||start) AS no_snp
    FROM variant_function group by file, region
    ORDER BY file, region;
")
snp_annot<- as.matrix(cbind(largew= snp_ann$no_snp[snp_ann$file == "var.flt.am_cat.vcf"],
                           lwl= snp_ann$no_snp[snp_ann$file == "var.flt.bmdm_cat.vcf"]))
rownames(snp_annot)<- unique(snp_ann$region)

arg.index<- match(c("exonic", "splicing", "ncRNA", "UTR5", "UTR3", "intronic", "upstream", "downstream", "intergenic"), rownames(snp_annot))

snp_annot<- snp_annot[order(snp_annot[,1]),]
windows(width=10/2.54, height= 8/2.54)
par(fg= 'grey25', cex= 0.85, mar= c(3,6,1,1), las=1, mgp= c(1.75, 0.75, 0), cex.main= 0.85)
barplot(t(snp_annot[arg.index,])/1000, beside= TRUE,
    col= c('dodgerblue', 'firebrick4'), horiz= TRUE, xlab= 'No. SNPs (x1000)',
    main= 'No. of SNPs by region')
legend('bottomright', bty= 'n', inset= c(0, -0.02), legend= c('L.w. x Land.', 'L. white'), col= c('firebrick4', 'dodgerblue'), pch= 19 )
savePlot('M:/Documents/LabBook/LabBook_Figures/20110415_snp_region.emf', 'emf')
graphics.off()


# ------------------------------------------------------------------------------
# exonic SNVs
# ------------------------------------------------------------------------------

exon.snv<- sqlQuery(conn, "
    select file, consequence, count(*) as no_snp
    from exonic_variant_function 
    group by file, consequence
    order by file, consequence;
")
exon.mut<- reshape(exon.snv, v.names= "no_snp", idvar= "consequence", timevar= "file", direction= 'wide')
colnames(exon.mut)<- c('consequence', 'largew', 'lwl')
## Make sure these names are correct!
exon.mut$consequence<- c('Non-synonymous', 'Stop gain', 'Stop loss', 'Synonymous', 'Frameshift del.', 'Frameshift ins.', 'Non-frameshift del.', 'Non-frameshift ins.')

## SNV in exonic regions
windows(width=12/2.54, height= 8/2.54)
nf <- layout( matrix(c(1,1,1,1,1,2,2,2) ))
par(fg= 'grey25', cex.axis=0.85, cex= 0.85, mar= c(1,8,2,1), las=1, mgp= c(1.75, 0.5, 0), cex.main= 0.85)
    barplot(t(exon.mut[which(exon.mut$lwl < 100),2:3]), beside= TRUE, xlim= c(0,60),
        col= c('dodgerblue', 'firebrick4'), horiz= TRUE, xlab= '',
        main= '', names.arg= exon.mut$consequence[which(exon.mut$lwl < 100)]
    )
    title(main= 'SNV in exonic regions')
    legend('topright', bty= 'n', inset= c(0, -0.02), legend= c('L.w. x Land.', 'L. white'), col= c('firebrick4', 'dodgerblue'), pch= 19 )
    par(mar= c(3,8,1.5,1))
    barplot(t(exon.mut[which(exon.mut$lwl > 100),2:3]), beside= TRUE,
        col= c('dodgerblue', 'firebrick4'), horiz= TRUE, xlab= 'No. SNVs',
        main= '', names.arg= exon.mut$consequence[which(exon.mut$lwl > 100)], xlim=c(0,8000)
    )
savePlot('M:/Documents/LabBook/LabBook_Figures/20110416_snp_exons.emf', 'emf')
graphics.off()

# ------------------------------------------------------------------------------
# Which genes/transcripts have variants?
# ------------------------------------------------------------------------------

vcf_ann<- sqlQuery(conn, 'select * from exonic_variant_function', stringsAsFactor= FALSE)

vcf_ann[1:10,]

vcf_annot<- strsplit(vcf_ann$annotation, ':')
vcf_annot2<- unlist(lapply(vcf_annot, strsplit, ',') )
vcf_transcripts<- unique(names(unlist(sapply(vcf_annot2, grep, patter= 'ENSSSCT',  perl= TRUE))))
vcf_transcripts[1:10]

length(vcf_transcripts)

## Annotate transcritps w/ biomaRt
library(biomaRt)
mart<- useDataset("sscrofa_gene_ensembl", useMart("ensembl"))
writeClipboard(listAttributes(mart)[1][,])
transcripts_annot<- getBM(	
             attributes= c("ensembl_transcript_id", "ensembl_gene_id", "external_transcript_id"), 
             filters= "ensembl_transcript_id",
             value= vcf_transcripts,
             mart= mart)
transcripts_annot$external_transcript_id[transcripts_annot$external_transcript_id == '']<- NA
sqlSave(conn, transcripts_annot, tablename= 'exonic_variant_function_genes', rownames= FALSE)
sqlQuery(conn, "ALTER TABLE exonic_variant_function_genes SET SCHEMA vcf_variants")
sqlQuery(conn, "COMMENT ON TABLE exonic_variant_function_genes IS 'Transcript names for the genes in exonic_variant_function.annotation. See 20110415_vcf_stats.R' ")
unique(transcripts_annot$external_transcript_id)

vcf_transcripts

# ------------------------------------------------------------------------------
# TRITUME
# ------------------------------------------------------------------------------

library(rpileup)
gcoord<- gene.coord(gene_id= 'ENSSSCG00000007585', genome= 'sscrofa') ## ACTB [ Ensembl gene: ENSSSCG00000007585 ]  

ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
mart<- useDataset("sscrofa_gene_ensembl", ensembl)

exon_annot<- getBM(	
             attributes= c("ensembl_transcript_id", "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end"), 
             filters= "ensembl_gene_id",
             value= 'ENSSSCG00000007585',
             mart= mart)

transcripts_annot<- getBM(	
             attributes= c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", "transcript_start", "transcript_end"), 
             filters= "ensembl_gene_id",
             value= 'ENSSSCG00000001404',
             mart= mart)


accepted_hits<- c("F:/data/bamfiles/20110202_rnaseq_am_ctrl/accepted_hits.bam",
                  "F:/data/bamfiles/20110202_rnaseq_am_lps/accepted_hits.bam",
                  "F:/data/bamfiles/20110202_rnaseq_bmdm_ctrl/accepted_hits.bam",
                  "F:/data/bamfiles/20110202_rnaseq_bmdm_lps/accepted_hits.bam")
accepted_hits_am<- c("F:/data/bamfiles/20110202_rnaseq_am_ctrl/accepted_hits.bam",
                  "F:/data/bamfiles/20110202_rnaseq_am_lps/accepted_hits.bam")


"7";29177013;29177013

ipile<- import.pileup(accepted_hits_am, rname='7', 29177013, 29177013, select= 'all' )
ipile<- import.pileup(accepted_hits_am, rname=gcoord[,1], gcoord[,2], gcoord[,3] )

abline(v= exon_annot$exon_chrom_start, col= 'grey', lty= 'dotted')
abline(v= exon_annot$exon_chrom_end, col= 'red', lty= 'dotted')

plot.pileup(ipile, overplot= FALSE, lwd= 2, type= 'l',
    col= c('dodgerblue', 'firebrick4'), pileup.names=c('am - ctrl', 'am - lps'))



plot.pileup(ipile, overplot= FALSE, lwd= 2, type= 'l',
    col= c('dodgerblue', 'firebrick4', 'dodgerblue', 'firebrick4'), pileup.names=c('am - ctrl', 'am - lps', 'bmdm - ctrl', 'bmdm - lps'))
for(i in 1:length(accepted_hits)){
    par(mfg=c(i,1))
    abline(v= exon_annot$exon_chrom_start, col= 'grey', lty= 'dotted')
    abline(v= exon_annot$exon_chrom_end, col= 'grey', lty= 'dotted')
}
par(mfg=c(2,1)); abline(v= exon_annot$exon_chrom_start)

transcripts_annot<- getBM(	
             attributes= c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", "transcript_start", "transcript_end"), 
             filters= "ensembl_gene_id",
             value= 'ENSSSCG00000001455',
             mart= mart)


"ENSSSCG00000001455"