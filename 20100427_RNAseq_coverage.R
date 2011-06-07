library(RODBC)
conn_db<- odbcConnect(dsn= 'pgCAGE', case='nochange', uid= 'dberaldi', pwd= 'mypassword')

# ---------------------------[ Using annotated GTF 28/4/10]-------------------------------

## This table created by 20100429_RNAseq_coverage.sql
transcripts<- sqlFetch(conn_db, 'tmp_transcript_cov', stringsAsFactors= FALSE)

## Get exon boundaries
exons<- sqlQuery(conn_db, stringsAsFactors= FALSE,
    "select distinct 
         array_to_string(regexp_matches(attributes, '.*transcript_id (.*?);.*'), '')::varchar AS transcript_id,
         array_to_string(regexp_matches(attributes, '.*gene_name (.*?);.*'), '|')::varchar AS gene_name,
         strand, f_start, f_end,
         array_to_string(regexp_matches(attributes, '.*exon_number (.*?);.*'), '|')::int AS exon_number
    from sus_scrofa_sscrofa9_56_gtf
    WHERE feature like 'exon'
    ORDER BY transcript_id, exon_number
    ")
odbcClose(conn_db)

transcripts[1:100,1:11]
names(transcripts)

(nctrl<- length(which(is.na(transcripts$ctrl_rpkm)==FALSE)))
(nlps<- length(which(is.na(transcripts$lps_rpkm)==FALSE)))

## windows(width= 20/2.54, height= 9/2.54)
par(mfrow= c(1,2), cex=0.75)
    ## Histogram of RPKMs
    hist(log2(transcripts$ctrl_rpkm), main= "Distribution of RPKM (ensembl transcripts)", xlab= 'Log(RPKM)', border= 'blue')
    hist(log2(transcripts$lps_rpkm), xlab= 'Log(RPKM)', border= 'red', add=T)
    abline(v= c(log2(2), log2(5)), col= "grey", lty= "dashed")
    legend('topleft', legend= c(paste('CTRL (n= ', nctrl, ')', sep=''), 
                                 paste('LPS (n= ', nlps, ')', sep=''), 'FPKM 2 and 5'
                                 ),
      col= c('blue', 'red', 'grey'), lty= 'solid', lwd= 2, y.intersp = 1.2, bty= 'n')
    
    ## --------------------[ Count of transcripts with rpkm >= x ]-------------------

    ntrans<- vector(length=0, mode= 'numeric')
    rpkm<- seq(0, 100)
    for(i in rpkm){
        ntrans<- append(ntrans, length(na.omit(transcripts$ctrl_rpkm[transcripts$ctrl_rpkm > i])))
        }
    plot(rpkm, ntrans, ylim= c(0, max(ntrans)), type= 'l', main= 'Number of transcripts/RPKM (CTRL)', xlab= 'RPKM', ylab= 'N. transcripts')
    points(x= c(0,2,5,10), y= ntrans[which(rpkm %in% c(0,2,5,10))], pch= 16, col= 'red')
    text(x= c(0,2,5,10) + 5, y= ntrans[which(rpkm %in% c(0,2,5,10))], 
        labels= paste('rpkm > ', c(0,2,5,10), '\nn= ', ntrans[which(rpkm %in% c(0,2,5,10))], sep=''), 
         adj= 0, cex= 1)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100427_RNAseq_coverage_Fig_1.emf', 'emf')



# -------------------------[ Read coverage for some genes  ]------------------------------

## set control / differentally expressed genes
gene<- c('ACTB_PIG', 'Q00P27_PIG', 'HPRT_PIG', 'TNFA_PIG', 'IL1A_PIG', 'IL8_PIG')



windows(width= 18/2.54, height= 20/2.54)
par(mfcol= c(3,2), cex=0.75, mar=c(0.25,3.5,0.25,2), oma= c(0,0,2,0))
for(g in gene){
      ctrl_t<- as.numeric(unlist(strsplit(transcripts[which(transcripts$gene_name == g), "ctrl_transcript_cov_5p3p"], ',')))
      lps_t<- as.numeric(unlist(strsplit(transcripts[which(transcripts$gene_name == g), "lps_transcript_cov_5p3p"], ',')))
      plot(x=1:length(ctrl_t), y= ctrl_t, type= 'l', ylim= c(0, max(c(ctrl_t, lps_t))), col= 'blue', xaxt= 'n', xlab='', 
           ylab='', bty= 'n')
      mtext(side=4, text= sub('_PIG', '', g))
      abline(h= 0)
      ex_size<- cumsum(exons[which(exons$gene_name == g), "f_end"] - exons[which(exons$gene_name == g), "f_start"])
      segments(x0= ex_size, y0= 0, y1= max(c(ctrl_t, lps_t)), lty='dashed', col= 'darkgrey', lwd= 1.5)
      points(x=1:length(lps_t), y= lps_t, type= 'l', col= 'red')
      arrows(x0= (length(lps_t) - 100), x1= length(lps_t), y0= max(c(ctrl_t, lps_t)) - max(c(ctrl_t, lps_t))*0.05, angle=90, code= 3, length= 0.05)
      text(x= (length(lps_t) - 120), y= max(c(ctrl_t, lps_t)) - max(c(ctrl_t, lps_t))*0.05, labels= '100 bp', adj= 1)
      }
title(main= 'Depth of coverage', outer= T)
mtext(side=2, text= 'n. reads', line=2.5, cex= 0.75)
legend(legend= c('CTRL', 'LPS'), 
    x= length(lps_t) - length(lps_t)*0.35, 
    y= (max(c(ctrl_t, lps_t)) - max(c(ctrl_t, lps_t))*0.05), 
    bty= 'n', lty= 'solid', col=c('blue', 'red'), yjust=0.5)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100427_RNAseq_coverage_Fig_2.emf', 'emf')
graphics.off()

# --------------------------------[ 3 and 5 prime coverage ]---------------------

## Extract first 100 bp (5') from each transcript

f_100<- function(xv){
    head5p<- 200
    tail3p<- 100
    ## Extract the 5' and 3' ends from a transcript represented as a string comma separated
    top<- head(as.numeric(strsplit(xv, ',')[[1]]), head5p)    ## 5' end
    top<- c(top, rep(NA, 200 - length(top)))               ## Pad the transcript with NA if shorter than the slice required
    bottom<- tail(as.numeric(strsplit(xv, ',')[[1]]), tail3p) ## 3' end
    bottom<- c(rep(NA, 100 - length(bottom)), bottom)
    return(c(top, bottom))
    }

ctrl_5p<- lapply(transcripts$ctrl_transcript_cov_5p3p, f_100)
names(ctrl_5p)<- transcripts$transcript_id
ctrl_5p<- t(do.call(cbind, ctrl_5p))
colnames(ctrl_5p)<- paste('ctrl_', 1:ncol(ctrl_5p), sep='')

lps_5p<- lapply(transcripts$lps_transcript_cov_5p3p, f_100)
names(lps_5p)<- transcripts$transcript_id
lps_5p<- t(do.call(cbind, lps_5p))
colnames(lps_5p)<- paste('lps_', 1:ncol(lps_5p), sep='')

names(transcripts)
lps_5rpkm<- lps_5p/(transcripts$lps_sum/transcripts$transcript_lenght)
std<- apply(lps_5rpkm, 2, sd, na.rm=T)
avg<- apply(lps_5rpkm, 2, mean, na.rm=T)
par(cex=0.95, mar=c(1,4,3,1), cex.main=1.1)
plot(1:100, avg[1:100], type= 'l', col= 'blue', xlim= c(1, 300), 
    ylim= c(0, max(avg)), lwd= 2, xlab= '', ylab= 'Relative coverage', xaxt='n', main= 'Average coverage vs position')
points(101:200, avg[101:200], type= 'l', col= 'green', lwd= 2)
points(201:300, avg[201:300], type= 'l', col= 'red', lwd= 2)
abline(v=c(100, 200), lty='dotted', col= 'darkgrey')
text(x= c(50, 150, 250), y= 0.7, labels=c("5' end\n(1:100 bp)", "5' ==> 3'\n(101:200 bp)", "3' end\n(last 100 bp)"), cex= 1)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100427_RNAseq_coverage_Fig_3.emf', 'emf')


arrows(x0=1:300, y0= avg-std, y1= avg+std, angle=90, length=0)
# boxplot(lps_5rpkm[sample(1:nrow(lps_5rpkm), 100),])
# boxplot(lps_5rpkm[1:200,])

lps_5rpkm[1:10,1:10]

lps_5rpkm<- subset(lps_5p, subset= (transcripts$lps_rpkm > 10) & (transcripts$lps_rpkm < 100), select= 1:100) ## & !is.na(transcripts$gene_name)
dim(lps_5rpkm)
set.seed(1234)
five<- lps_5rpkm[sample(1:nrow(lps_5rpkm), 20),]
boxplot(five)
plot(x= 1:100, y= five[1,], type='n', ylim=c(0, max(five, na.rm=T)))
for(i in 1:nrow(five)){
    d<- five[i,]    
    points(x= 1:100, y= d, col= i, type= 'l')
    }

## test if there is a decline when nearing the 5' end
lps_res<- lps_5p[transcripts$lps_rpkm > 10, ]
tr<- rep(rownames(lps_res), each= ncol(lps_res))
pos<- rep(colnames(lps_res), nrow(lps_res))
cover<- as.vector(t(lps_res))
lps_lm<- data.frame(tr, pos, cover)

lm1<- lm(cover ~ pos + tr, data= lps_lm)
f_100(transcripts$lps_transcript_cov_5p3p[1])

dim(lps_lm)

as.numeric(strsplit(transcripts$ctrl_transcript_cov_5p3p[1], ',')[[1]][1:100])


## ------------------------------[ Tritume ]-------------------------------------
## Select transcript with RPKM > 5 in either lib
exons_rpkm5<- subset(exons, 
    ((is.na(exons$ctrl_rpkm) == FALSE & exons$ctrl_rpkm > 5) | 
    (is.na(exons$lps_rpkm) == FALSE & exons$lps_rpkm > 5)) & 
     exons$exon_number == 1, 
     select= 1:ncol(exons))
dim(exons_rpkm5)
exons_rpkm5[1:10, ]

yv<- sapply(as.character(exons_rpkm5$ctrl_cov_5p3p), function(x) as.numeric(unlist(strsplit(x, ','))))
yv<- lapply(yv, function(x) c(x, rep(NA, 100 - length(x))))
names(yv)<- exons_rpkm5$transcript_id
yv<- t(do.call(cbind, yv))
colnames(yv)<- paste('ctrl_', 1:ncol(yv), sep='')

xv<- sapply(as.character(exons_rpkm5$lps_cov_5p3p), function(x) as.numeric(unlist(strsplit(x, ','))))
xv<- lapply(xv, function(x) c(x, rep(NA, 100 - length(x))))
names(xv)<- exons_rpkm5$transcript_id
xv<- t(do.call(cbind, xv))
colnames(xv)<- paste('lps_', 1:ncol(xv), sep='')

exons_rpkm5<- cbind(exons_rpkm5[1:10], yv, xv)

plot(1:100, exons_rpkm5[5, 11:110], type= 'l')

boxplot(exons_rpkm5[, 10:ncol(exons_rpkm5)])

exons_rpkm5[5,1:10 ]

# -------------------------------------------------------------------------------------------


ref_trans<- sqlFetch(conn_db, 'tmp_ref_trans')
ref_trans
exons<- sqlFetch(conn_db, 'tmp_exons')
head(exons)



par(mfrow= c(3,2), mar= c(0.5,4,0.5,1), cex=0.8, oma= c(0,0,2,0))
  for(trans_id in unique(exons$transcript_id)){
    with(exons[exons$transcript_id== trans_id,], {
        plot(exon_offset, coverage, type= 'l', lwd= 2, ylab= 'Coverage (# reads)', xaxt='n', xlab= '', frame.plot=F)
        abline(h= 6, col= 'blue', lty='dotted')
        segments(x0= cumsum(ref_trans$feature_length[ref_trans$transcript_id == trans_id]), y0= 0, y1= max(coverage), col= 'red', lty= 'dotted')
        mtext(side= 3, line= -1, text= gene_name[1], cex= 1, adj= 0.85, font=3)
        text(x= 0, y= max(coverage)- max(coverage)*0.05, labels= '100 bp', adj= 0)
        arrows(x0= 0, x1= 100, y0= max(coverage)- max(coverage)*0.1, lwd= 2, angle= 90, length= 0.05, code=3)
        })
    }
  with(exons[exons$transcript_id== trans_id,], {
      })
title(main= 'Read coverage for some genes of interest', outer= TRUE)

# ---------------------------[ Using annotated GTF 28/4/10]-------------------------------

exons<- sqlFetch(conn_db, 'tmp_ss9_56_gtf_coverage')
dim(exons[exons$lps_nreads != '' & exons$ctrl_nreads != '', ])

yv<- sapply(exons$ctrl_nreads, function(x) as.numeric(unlist(strsplit(as.character(x), ','))))
yv<- lapply(yv, function(x) c(x, rep(NA, 100 - length(x))))
yv<- t(do.call(cbind, yv))
colnames(yv)<- paste('ctrl_', 1:ncol(yv), sep='')

xv<- sapply(exons$lps_nreads, function(x) as.numeric(unlist(strsplit(as.character(x), ','))))
xv<- lapply(xv, function(x) c(x, rep(NA, 100 - length(x))))
xv<- t(do.call(cbind, xv))
colnames(xv)<- paste('lps_', 1:ncol(xv), sep='')


# -------------------------------------[ Tritume ]------------------------------------

exons<- sqlQuery(conn_db, 'SELECT * FROM tmp_ss9_56_gtf_coverage ORDER BY transcript_id, exon_number')
head(exons)
exons_rpkm<- subset(exons, !is.na(exons$ctrl_rpkm) | !is.na(exons$lps_rpkm), select= 1:ncol(exons))

f_exon<- exons_rpkm[exons_rpkm$exon_number == 1, ]
dim(f_exon)



f_exon<- cbind(f_exon[, 1:5], yv, xv)

f_exon[1:10,1:10]

dim(f_exon)
dim(xv)

lps_nreads<- xv

dim(y)
f_exon[1:10,]
yv[1:10,1:10]
dim(yv)

lapply(xv[1:10], length)
length(xv[1])


xv
dim(f_exon)


xv[1:10]
strsplit(as.character(f_exon$lps_nreads[1]), ',')
f_exon$lps_nreads[1:10]

head(f_exon)

(m<- t(matrix(c(10,10,10,10,2,15,15,15,15,5,24,24,24,24,12), ncol=3)))
(d<- c(2, 5, 12))
(m/d)

