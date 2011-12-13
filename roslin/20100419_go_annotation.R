 _____________________________________________________________________________
#
# GO analysis of ensembl transcripts from RNAseq.
# Determine whether transcripts differemtially expressed between CTRL and LPS
# have some GO terms over-represented.
#
# Using library topGO
#
# Data: a table with transcript IDs and some measure of significance of
# differential expression (p-value, logFC)
# 
# GO database: A table with one transcript/row and the associated GO terms 
#
# ____________________________________________________________________________

library(topGO)

# -------------------------[ User's input ]-----------------------------------

## Transcripts from RNAseq
gene_universe<- read.table('U:/Documents/GO_annotation/rnaseq_gene_universe_fpkm_5.txt', sep= '\t', header=TRUE)

## GO graph - GO annotation for all EMSEMBL transcripts, S. scrofa 57 (transcripts w/o any GO removed)
ensembl2GO<- readMappings('U:/Documents/GO_annotation/go_annotation_ensembl_sscrofa57.txt', sep= '\t', IDsep= ',')

# -----------------------[ Prepare topGOdata ]--------------------------------

## Define *gene universe*: 
## All the transcripts in the libraries filtered to have:
## - Some GO term associated
## - A count in both LPS and CTRL
## - logFC >x (or some other metrics to detect some variation between treatments)
##
dim(gene_universe)

## From the gene universe create a named numeric vector.
## The numeric value is the *score*, i.e. the measure of significance for diff. expr. (e.g. p-value)
## The name is the transcript ID
##
gene_list<- gene_universe$FDR
names(gene_list)<- gene_universe$transcript_id

## Define a function that selects from the gene_list (universe) the transcripts for which
## GO over-representation has to be tested. Typically, the genes with high sig. p-value or large logFC
## 
select_genes<- function(gene_list){
    ## Returns TRUE/FALSE for the genes satisfying the conditions
    return(gene_list < 0.01 )    ## logFC select genes according to the direction of effects & gene_universe$logFC > 0
    } 

## Some summary for the gene universe and selected list
## 
n_geneu<- length(gene_list)
n_int<- length(which(select_genes(gene_list) == TRUE)) ## Count of interesting Genes returned by filtering
sel_perc<- format(n_int/length(gene_list)*100, digits= 2)
cbind('Gene universe (n)'= n_geneu, 'Tested genes (n)'= paste(n_int, ' (', sel_perc, '%)',  sep='' ))

## Compile topGOdata object
##
GOdata<- new('topGOdata', # object class name
  nodeSize= 10,            # Min number of genes associated which each GO term
  ontology= 'BP',         # GO category (MF, BP, CC)
  allGenes= gene_list,    # Named, numeric vector representing universe
  geneSel= select_genes,  # The function that defines which genes in allGenes are 'interesting'
  annot= annFUN.gene2GO,  # Function to link GO terms in ensembl2GO to transcript IDs
  gene2GO= ensembl2GO     # Table of transcripts IDs and their GO terms
  )
GOdata

# ------------------------------[ Analyse data ]-----------------------------------

## Perform tests
## 
(resultKS<- runTest(GOdata, algorithm= 'elim', statistic= "ks"))

# ------------------------------[ Visualize results ]------------------------------

use_test<- resultKS
geneData(use_test)
(allRes<- GenTable(GOdata, KS= use_test, topNodes= 50)) ## View top scoring GO terms
write.table(allRes, file= 'C:/Tritume/GOdata_bp.txt', sep= '\t', row.names=F)

## GO graph
setwd('M:/Documents/LabBook/LabBook_Figures')

showSigOfNodes(GOdata, score(use_test), firstSigNodes = 5, useInfo = 'def') ##, wantedNodes= 'GO:0045649'
savePlot(file= "20100419_go_annotation.R.Fig_1.emf", type='emf')
dev.off()
printGraph(GOdata, use_test, firstSigNodes = 5, fn.prefix = "20100419_go_annotation.R.Fig_1", useInfo = "def", pdfSW = TRUE)

showGroupDensity(GOdata, whichGO= "GO:0005125")

## Genes in specific GO term(s)
sel.terms<- 'GO:0006955'
ann.genes <- genesInTerm(GOdata, sel.terms)
str(ann.genes)
length(unlist(ann.genes))

ann.score <- scoresInTerm(GOdata, sel.terms)
str(ann.score)
length(unlist(ann.score))

genes_score<- data.frame(ann.genes= unlist(ann.genes), ann.score= unlist(ann.score))
genes_score[1:10, ]
dim(genes_score[genes_score$ann.score < 0.01,])


num.ann.genes <- countGenesInTerm(GOdata)
str(num.ann.genes)



termStat(GOdata, sel.terms)


# -----------------------------------[ Tritume ]-----------------------------------
myInterstedGenes<- gene_universe$transcript_id[gene_universe$PValue<0.01]
geneList<- factor(as.integer(gene_universe$transcript_id %in% myInterstedGenes))
names(geneList)<- gene_universe$transcript_id
str(geneList)

#library(ALL)
#data(ALL)
#data(geneList)
#hist(geneList)

x<- 20                   ## number of white balls (interesting genes) drawn without replacement from an urn which contains both black and white balls
m<- 126                 ## the number of white balls (interesting genes) in the urn
n<- 4193 - m            ## the number of black balls (gene universe - interestin genes) in the urn.
k<- 533                 ## the number of balls drawn from the urn (genes annotated to each GO term).
dhyper(x, m, n, k)

5              (126 - 5)= 121
(102 - 5)=97  (4193 - (5 + 121 + 97))=3970

x<- 20     ## Number of interesting genes drawn at given GO term (white balls drawn)
m<- 126    ## Total number of intersting genes (number of white balls in the urn)
k<- 533    ## Number of genes (interesting and non) annotated at given GO term (num. of balls drawn from urn)
n<- 4193   ## Total number number of genes in the universe (black + white balls in the urn)

## Fisher's exact test:
x                        ## Same as above
wb<- m - x               ## White balls remaining in the urn after the draw 
bb<- k - x               ## Black balls in the draw
nn<- n - sum(c(x,wb,bb))

data_xm<- c(x, wb, bb, nn)
sum(data_xm)  ## Same as n

(xm<- matrix(data= data_xm, nrow= 2, byrow= T))
fisher.test(xm)

goID <- "GO:0006091"


go <- usedGO(GOdata) ## the available GO terms (all the nodes in the graph)

rnaseq_data[1:10,]

dim(gene_universe)
