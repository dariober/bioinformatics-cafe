#
#  Extract the first n bases of the TSS from the human genome
#  for the CAGE tags in the FANTOM5 library Macrophage%20-%20monocyte%20derived%20donor1.CNhs10861.11232-116C8.hg19.ctss.bed.gz
#
#  Use BSgenome utility
#

library(RODBC)
library(BSgenome)

# Download Hsapiens genome (do it only once)
#
## source("http://bioconductor.org/biocLite.R")
## biocLite("BSgenome.Hsapiens.UCSC.hg19")

library(BSgenome.Hsapiens.UCSC.hg19)

## Get the first 30 bases of the from bed file of the TSSs from the FANTOM5 macrophage library.
## BED file is in 0-based coordinates so the tss is on f_start + 1.
conn<- odbcConnect(dsn= 'pgVitelleschi')
mac_tss<- sqlQuery(conn, "
    SELECT ctss_tag, rname, CASE WHEN strand = '+' THEN (f_start + 1)
                                 WHEN strand = '-' THEN (f_start + 1) - 29 END AS f_start,
                            CASE WHEN strand = '+' THEN (f_start + 1) + 29
                                 WHEN strand = '-' THEN (f_start + 1)      END AS f_end,
           strand,
           tpm
    FROM macrophage_monocyte_derived_hg19_ctss_bed
    ORDER BY rname, f_start;
    ")
odbcClose(conn)
mac_tss[1:10,]
# Extract sequences
tss_seq<- getSeq(x= Hsapiens, names= mac_tss$rname, start= mac_tss$f_start, end= mac_tss$f_end, strand= mac_tss$strand)

# Bind sequence to tags adn write out
mac_tss<- cbind(mac_tss, tss_seq)

write.table(mac_tss, 'F:/data/20110105_FANTOM5/Macrophage%20-%20monocyte%20derived%20donor1.CNhs10861.11232-116C8.hg19.ctss.seq.txt', sep='\t', row.names=F, quote=F, col.names=T)

" -------------------------------[ TRITUME ]------------------------------- "

ctss_tag_id<- data.frame(serial= as.numeric(seq(1, length(tss_seq))), ctss_tag= paste('>', as.character(mac_tss$ctss_tag), sep=''))
ctss_seq<- data.frame(serial= as.numeric(seq(1, length(tss_seq))), ctss_seq= as.character(mac_tss$tss_seq)

fasta.1<- rbind(ctss_tag_id, ctss_seq)
fasta.1<- fasta.1[order(fasta.1$serial),]

tss_fasta<- 
object.size(tss_seq)/10^6

mac_tss[which(mac_tss$ctss_tag == 'chr7:5570231..5570232,-'), ]


> all.genomic<-getSeq(Hsapiens, the.chrom, starts, ends)
>
> where
> the.chrom[1:5]
> [1] "chr1" "chr1" "chr1" "chr1" "chr1"'
> starts[1:5]
> 3187526 3487463 3777276 4144186 4274111
> > ends[1:5]
> [1] 3187790 3487763 3777555 4144499 4274416
>
>
> etc etc...
> myb<-"YAACKG"
> length(all.genomic)
> system.time(x<- XStringViews(all.genomic, "DNAString"))
> x.labels<-paste(the.chrom,starts,ends,sep=":")
> names(x)<-x.labels
> ###################### forward counts #################
> all.matches<-matchPattern(myb.dna,x,max.mismatch=0, fixed=FALSE) # needs
> a stringView to vectorize
> the.cov<-coverage(all.matches)
>
> counts<-aggregate(the.cov,start=start(x),end=end(x),FUN=sum)/length(myb.dna)


conn<- odbcConnect(dsn='pgVitelleschi')
tnf<- sqlQuery(conn, "
    select * from macrophage_monocyte_derived_hg19_ctss_bed 
    where rname like 'chr6' and f_start between (31543344-1000) and (31543344+1000)
    ")
plot(x= tnf$f_start, y= tnf$tag_count, type= 'h')


library("BSgenome.Mmusculus.UCSC.mm9")



mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
getSequence(chromosome=6, start= 31543342-10, end= 31543343+10, seqType='5utr', mart=mart)