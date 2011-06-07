# -----------------------------------------------------------------------------
#  Fetch from biomaRt attributes for genes of interest
# -----------------------------------------------------------------------------

library(biomaRt)

## Connect to Biomart: Species to use
mart_hs<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mart_mm<- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
mart_ss<- useDataset("sscrofa_gene_ensembl", useMart("ensembl"))

# -----------------------------------------------------------------------------
# Genes of interest
# -----------------------------------------------------------------------------

# Add here more gene_ids or more species *AND* reference gene names. Always in the same order!!
gene_list_hsapiens<- c("ENSG00000102575", "ENSG00000042493", "ENSG00000115009", "ENSG00000172216", "ENSG00000182578", "ENSG00000148180", "ENSG00000131203", "ENSG00000017427", "ENSG00000168811", "ENSG00000113302", "ENSG00000125538", "ENSG00000136244", "ENSG00000187098", "ENSG00000066336", "ENSG00000138378", "ENSG00000232810", "ENSG00000075624")
gene_list_mmusculus<- c("ENSMUSG00000001348", "ENSMUSG00000056737", "ENSMUSG00000026166", "ENSMUSG00000056501", "ENSMUSG00000024621", "ENSMUSG00000026879", "ENSMUSG00000031551", "ENSMUSG00000020053", "ENSMUSG00000027776", "ENSMUSG00000004296", "ENSMUSG00000027398", "ENSMUSG00000025746", "ENSMUSG00000035158", "ENSMUSG00000002111", "ENSMUSG00000062939", "ENSMUSG00000024401", "ENSMUSG00000029580")
gene_list_sscrofa<- c("ENSSSCG00000013612", "ENSSSCG00000008239", "ENSSSCG00000016254", "ENSSSCG00000007468", "ENSSSCG00000014441", "ENSSSCG00000005515", "ENSSSCG00000007007", "ENSSSCG00000000857", "ENSSSCG00000011730", "ENSSSCG00000017044", "ENSSSCG00000008088", "ENSSSCG00000015385", "ENSSSCG00000011512", "ENSSSCG00000013237", "ENSSSCG00000016059", "ENSSSCG00000001404", "ENSSSCG00000007585")
reference_gene_names<- c("ACP5", "CAPG", "CCL20", "CEBPB", "CSF1R", "GSN", "IDO1", "IGF1", "IL12A", "IL12B", "IL1B", "IL6", "MITF", "SPI1", "STAT4", "TNF", "ACTB")

## This table binds the vectors above. 
ref_table<- data.frame(cbind(ensembl_gene_id= c(gene_list_hsapiens, gene_list_mmusculus, gene_list_sscrofa), reference_gene_name= rep(reference_gene_names, 3)))

## Fetch attributes: TRANSCRIPTS
trans_attr<- data.frame(rbind(
trans_hs<-getBM(
    mart= mart_hs, filters= c("ensembl_gene_id"), values= gene_list_hsapiens, 
    attributes= c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", "end_position", "strand",
                  "ensembl_transcript_id", "external_transcript_id", "transcript_start", "transcript_end"))
trans_hs$species<- "hsapiens"

trans_mm<- getBM(mart= mart_mm, 
    filters= c("ensembl_gene_id"), values= gene_list_mmusculus, 
    attributes= c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", "end_position", "strand",
                 "ensembl_transcript_id", "external_transcript_id", "transcript_start", "transcript_end"))
trans_mm$species<- 'mmusculus'
                 
trans_ss<- getBM(mart= mart_ss, 
    filters= c("ensembl_gene_id"), values= gene_list_sscrofa, 
        attributes= c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", "end_position", "strand",
                    "ensembl_transcript_id", "external_transcript_id", "transcript_start", "transcript_end"))
trans_ss$species<- 'sscrofa'

## Bind all together
trans_attr<- rbind(trans_hs, trans_mm, trans_ss)

## Add reference_gene_name column
trans_attr<- merge(trans_attr, ref_table, by= intersect("ensembl_gene_id", "ensembl_gene_id"))
## Replace strand 1/-1 with +/- for consistency with BED format
trans_attr$strand[trans_attr$strand == 1]<- '+'
trans_attr$strand[trans_attr$strand == -1]<- '-'

## This table uploaded to postgres immunogenes as transcript_attributes.
write.table(trans_attr, file= 'D:/Tritume/transcript_attributes.txt', sep='\t', row.names=F)

# -----------------------------------------------------------------------------
# The code from here below is probably deprecated
# -----------------------------------------------------------------------------

## Fetch attributes: GENES -- This block deprecated!
gene_attr<- data.frame(rbind(
    getBM(mart= mart_hs, 
      filters= c("ensembl_gene_id"),
      values= gene_list_hsapiens, 
      attributes= c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", "end_position", "strand")),
    getBM(mart= mart_mm, 
      filters= c("ensembl_gene_id"),
      values= gene_list_mmusculus, 
      attributes= c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", "end_position", "strand")),
    getBM(mart= mart_ss, 
      filters= c("ensembl_gene_id"),
      values= gene_list_sscrofa, 
      attributes= c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", "end_position", "strand"))
    ))
gene_attr$species<- rep(c('hsapiens', 'mmusculus', 'sscrofa'), each= c(length(gene_list_hsapiens), length(gene_list_mmusculus), length(gene_list_sscrofa)))
## Add reference_gene_name column
gene_attr<- merge(gene_attr, ref_table, by= intersect("ensembl_gene_id", "ensembl_gene_id"))
## Replace strand 1/-1 with +/- for consistency with BED format
gene_attr$strand[gene_attr$strand == 1]<- '+'
gene_attr$strand[gene_attr$strand == -1]<- '-'
write.table(gene_attr, file= 'D:/Tritume/gene_attributes.txt', sep='\t', row.names=F)

# -----------------------------------------------------------------------------
# Get the full gene sequences +300bp upstream 
# -----------------------------------------------------------------------------

## Humans only upstream seqs.
write.table(
  cbind (
        species= 'hsapiens', upstream= 300,
        getSequence(id= gene_list_hsapiens, type= "ensembl_gene_id", seqType= "gene_exon_intron", upstream= 300, mart= mart_hs)
  ), 
  file= 'D:/Tritume/gene_exon_intron_up300_hs.txt', sep='\t', row.names=F, quote= FALSE)

## Mouse
write.table(
  cbind (
        species= 'mmusculus', upstream= 300,
        getSequence(id= gene_list_mmusculus, type= "ensembl_gene_id", seqType= "gene_exon_intron", upstream= 300, mart= mart_mm)
  ), 
  file= 'D:/Tritume/gene_exon_intron_up300_mm.txt', sep='\t', row.names=F, quote= FALSE)

## Pig
write.table(
  cbind (
        species= 'sscrofa', upstream= 300,
        getSequence(id= gene_list_sscrofa, type= "ensembl_gene_id", seqType= "gene_exon_intron", upstream= 300, mart= mart_ss)
  ), 
  file= 'D:/Tritume/gene_exon_intron_up300_ss.txt', sep='\t', row.names=F, quote= FALSE)


# -----------------------------------------------------------------------------
# Try to get sequences flanking coding genes
# -----------------------------------------------------------------------------

## Humans only upstream seqs.
write.table( cbind(
    species= 'hsapiens', upstream= 300,
    getSequence(id= gene_list_hsapiens, type= "ensembl_gene_id", seqType= "coding_gene_flank", upstream= 300, mart= mart_hs)
  ),
  file= 'D:/Tritume/coding_gene_flank_up300_hs.txt', sep='\t', row.names=F, quote= F)

## Mouse
write.table( cbind(
    species= 'mmusculus', upstream= 300,
    getSequence(id= gene_list_mmusculus, type= "ensembl_gene_id", seqType= "coding_gene_flank", upstream= 300, mart= mart_mm)
  ),
  file= 'D:/Tritume/coding_gene_flank_up300_mm.txt', sep='\t', row.names=F, quote= F)

## Pig
write.table( cbind(
    species= 'sscrofa', upstream= 300,
    getSequence(id= gene_list_sscrofa, type= "ensembl_gene_id", seqType= "coding_gene_flank", upstream= 300, mart= mart_ss)
  ),
  file= 'D:/Tritume/coding_gene_flank_up300_ss.txt', sep='\t', row.names=F, quote= F)




# -----------------------------------------------------------------------------
# TRITUME
# -----------------------------------------------------------------------------sequences<- rbind(

reference_genes_norm<- cbind(
    reference_gene_name= rep(reference_gene_names, 3),
    ensembl_gene_id= c(gene_list_hsapiens, gene_list_mmusculus, gene_list_sscrofa),
    species= rep(c('hsapiens', 'mmusculus', 'sscrofa'), each= c(length(gene_list_hsapiens), length(gene_list_mmusculus), length(gene_list_sscrofa)))
    );
write.table(reference_genes_norm, file= 'D:/Tritume/reference_genes_norm.txt', sep='\t', row.names=F)

reference_genes<- data.frame(cbind(reference_gene_name= reference_gene_names, ensembl_gene_id_hs= gene_list_hsapiens, ensembl_gene_id_mm= gene_list_mmusculus, ensembl_gene_id_ss= gene_list_sscrofa))
gene_description<- getBM(mart= mart_hs, filters= c("ensembl_gene_id"), 
    attributes= c("ensembl_gene_id", "description"), values= gene_list_hsapiens)
names(gene_description)[1]<- "ensembl_gene_id_hs"
reference_genes_descr<- merge(reference_genes, gene_description, by= intersect("ensembl_gene_id_hs",  "ensembl_gene_id_hs"))
reference_genes_descr<- reference_genes_descr[,c("reference_gene_name", "ensembl_gene_id_hs",  "ensembl_gene_id_mm",  "ensembl_gene_id_ss", "description")]
write.table(reference_genes_descr, file= 'D:/Tritume/reference_genes_ct.txt', sep='\t', row.names=F)


write.table(
  getSequence(id= gene_list_hsapiens, type= "ensembl_gene_id", seqType= "coding_gene_flank", downstream= 100, mart= mart_hs),
  file= 'D:/Tritume/coding_gene_flank_down100_hs.txt', sep='\t', row.names=F)


write.table(
  getSequence(id= gene_list_sscrofa, type= "ensembl_gene_id", seqType= "coding_gene_flank", downstream= 100, mart= mart_ss),
  file= 'D:/Tritume/coding_gene_flank_down100_ss.txt', sep='\t', row.names=F)

write.table(
  getSequence(id= gene_list_mmusculus, type= "ensembl_gene_id", seqType= "coding_gene_flank", downstream= 100, mart= mart_mm),
  file= 'D:/Tritume/coding_gene_flank_down100_mm.txt', sep='\t', row.names=F)



    getSequence(mart= mart_hs, id= gene_list_hsapiens, type= "ensembl_gene_id", seqType= "gene_exon_intron"),
    getSequence(mart= mart_mm, id= gene_list_mmusculus, type= "ensembl_gene_id", seqType= "gene_exon_intron"),
    getSequence(mart= mart_ss, id= gene_list_sscrofa, type= "ensembl_gene_id", seqType= "gene_exon_intron")
    )
sequences$type<- 'gene_exon_intron'



sequences_upstream<- data.frame("coding_gene_flank"= NULL, "ensembl_gene_id"=NULL, "type"=NULL)
for(i in (1:length(gene_list_hsapiens))){
    print(i)
    seq<- getSequence(id= gene_list_hsapiens[2], type= "ensembl_gene_id", seqType= "coding_gene_flank", upstream= 300, mart= mart_hs)
    seq$type<- "hsapiens"
    sequences_upstream<- rbind(sequences_upstream, seq)
    }

sequences_upstream<- rbind(
    getSequence(id= gene_list_hsapiens, type= "ensembl_gene_id", seqType= "coding_gene_flank", upstream= 300, mart= mart_hs), 
    getSequence(mart= mart_mm, id= "ENSMUSG00000031551", type= "ensembl_gene_id", seqType= "coding_gene_flank", upstream= 300),
    getSequence(mart= mart_ss, id= gene_list_sscrofa, type= "ensembl_gene_id", seqType= "coding_gene_flank", upstream= 300)
    ) ## Note sometimes an error like below appears. Repeat the query a few times and it will work.
      ## Error in getBM(c(seqType, type), filters = c(type, "upstream_flank"),  : 
      ## Query ERROR: caught BioMart::Exception::Usage: Filter upstream_flank NOT FOUND
sequences_upstream$type<- "coding_gene_flank"
write.table(sequences_upstream, file= 'D:/Tritume/sequences_upstream.txt', sep='\t', row.names=F)
aggregate(sequences_upstream$ensembl_gene_id, list(sequences_upstream$ensembl_gene_id), length)


entrez = c("673", "7157", "837")
getSequence(id = entrez, type = "entrezgene", seqType = "coding_gene_flank", upstream = 100, mart = mart_hs)
getSequence(id = gene_list_hsapiens, type = "ensembl_gene_id", seqType = "coding_gene_flank", upstream = 300, mart = mart_hs)


write.table(sequences, file= 'D:/Tritume/sequences.txt', sep='\t', row.names=F)

sequence_hs[1:5,]

library(biomaRt)
mart_hs<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# gene_list_hsapiens<- "..."
ref_genes<- getBM(mart= mart_hs,
    filters= c("ensembl_gene_id"),
    values= gene_list_hsapiens, 
    attributes= c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", "end_position", "strand", "description"))
write.table(ref_genes, 'D:/Tritume/ref_genes.txt', row.names= F, sep= '\t', quote=F)
  
  
exportFASTA(
    getSequence(id= "ENST00000373818", type= "ensembl_transcript_id", seqType= "coding", mart= mart_hs),
    'D:/Tritume/coding_gsn.fa')

exportFASTA(
    getSequence(id= "ENSSSCT00000006066", type= "ensembl_transcript_id", seqType= "coding", mart= mart_ss),
    'D:/Tritume/coding_gsn.fa')
  
exportFASTA(
    getSequence(id= "ENSMUST00000028239", type= "ensembl_transcript_id", seqType= "coding", mart= mart_mm),
    'D:/Tritume/coding_gsn.fa')
  