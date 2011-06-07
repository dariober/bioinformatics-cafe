#_________________________________________________________________________
#
#  GO analysis of RNAseq LPS and CTRL libraries. Trannscripts assembled 
#  by cufflinks aided by GTF file.
#
#  Using package GOstats
#_____________________________________________________________________________

library(GOstats)
library(GSEABase)
library(Rgraphviz)
library(gplots)

setwd('U:/Documents/GO_annotation/GOstats/output')

# -----------------------------[ Prepare GO database ]------------------------

## Load GO annotation prepared with 20100422_go_annotation_GOstats.sql
goframeData<- read.table('U:/Documents/GO_annotation/GOstats/go_annotation_ensembl_sscrofa57_GOstats.txt', sep='\t', header= TRUE)
goframeData$evidence<- as.character(goframeData$evidence)

goFrame<- GOFrame(goframeData, organism = "Sus scrofa")
goAllFrame = GOAllFrame(goFrame)

gsc<- GeneSetCollection(goAllFrame, setType = GOCollection())

## GO database from library(GO.db)
go_db<- as.list(GOTERM)

# --------------------[ Define gene universe and genes of interest]------------

## Transcript from RMNAseq with having GO annotation and 
## combined FPKM (ctrl + lps) > 5 (exclude low expressed genes)
##
gene_universe<- read.table('U:/Documents/GO_annotation/rnaseq_gene_universe_fpkm_5.txt', sep= '\t', header=TRUE)
genes<- gene_universe$transcript_id[gene_universe$FDR < 0.01]  ## Genes of interest (white balls); param: geneIds
universe<- gene_universe$transcript_id

#______________________________________________________________________________
#
#                         Analyse ontologies
#
#                         MOLECULAR FUNCTION
#______________________________________________________________________________


## Compile GOstats object
##
params <- GSEAGOHyperGParams(name = "S. scrofa custom GSEA based annot Params",
  geneSetCollection = gsc, geneIds = genes, universeGeneIds = universe,
  ontology = "MF", pvalueCutoff = 0.05, conditional = FALSE,
  testDirection = "over")

## Perform test(s)
## 
MF_over_condF <- hyperGTest(params)

# -------------------------------[ MF Visualization ]--------------------------

(test_summary<- summary(MF_over_condF, pvalue= 0.01, categorySize= 10))

htmlReport(MF_over_condF, file= 'mf_over_condF.html')

##
##  Plot GO terms
## 

top_n<- c(1:10) ## Lines in test_summary to plot in the graph
topGo<- test_summary$GOMFID[top_n] ## GO selected
gog<- GOGraph(topGo, GOMFPARENTS)  ## Prepare graph object
nodes_go<- nodes(gog) ## GO IDs in graph gog

## Convert GO accession to GO name
nodes_name<- vector(mode='character')
for(n in nodes_go){
    go_acc<- go_db[n]
    nodes_name<- append(nodes_name, Term(go_db[[n]]))
    }

## Make nodes redder if more significant (besed on pvalue: GO terms should be ordered by pvalue maybe check this)
##
alpha_gradient<- heat.colors(n=2, rep(seq(1, 0.3, length.out= 4), each= ceiling(length(top_n)/4)))[c(seq(1,10,by=2), seq(2,10, by=2))]
fillcol<- rep('lightblue', length(nodes_go)) ## Colour for GO terms not in the selected list (i.e. parents of thje selected terms)
fillcol[nodes_go %in% topGo]<- alpha_gradient

#plot(1:10, pch=16, cex= 5, col= cc[c(seq(1,10,by=2), seq(2,10, by=2))]) #<-- Test color palette

## see ?GraphvizAttributes for more attributes
node.attrs <-
    makeNodeAttrs(gog, shape="ellipse",
                  label = substr(nodes_name,1,12), ## Use substr() to shorten the names
                  fillcolor = fillcol,
                  fontsize= 24.0, ## Has no effect!
                  width= 2      ## Width of the nodes

                  )
windows(width= 17/2.54, height= 17/2.54)
plot(agopen(gog,
            recipEdges="distinct", layoutType="dot",
            nodeAttrs = node.attrs, name = "foo"))
# locator()
text(x= 0, y= 230, labels= "GO: Molecular function\nTest for over-representation", adj=0)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100422_go_annotation_GOstats_Fig_1.jpeg', 'jpeg')

# ---------------------------[ MF Genes in GO terms ]--------------------------

head(gene_universe)
topGenes<- unlist(geneIdsByCategory(MF_over_condF)[topGo[1:10]])  ## Returns the 'interesting' gene IDs on each GO term
go_mf_topgenes<- gene_universe[gene_universe$transcript_id %in% topGenes & gene_universe$associated_gene_name != '', ]
length(unique(go_mf_topgenes$associated_gene_name))


locus1<- vector(mode= 'character', length=0)
locus2<- vector(mode= 'character', length=0)
value<- vector(mode= 'numeric', length=0)

g1<- as.character(go_mf_topgenes$associated_gene_name)
logfc<- abs(go_mf_topgenes$logFC)
for(i in seq(1, nrow(go_mf_topgenes)) ){
    locus1<- append(locus1, g1)
    locus2<- append(locus2, rep(g1[i], length(g1)))
    value<- append(value,  1/abs((logfc - logfc[i])))
    }
exprmat<- data.frame(locus1, locus2) ## value
exprmat<- exprmat[exprmat$locus1 != exprmat$locus2, ]

write.table(exprmat, 'C:/Tritume/go_biolayout.txt', sep='\t', row.names=F)

## Pretty table of differentially expressed genes (topGenes) in 'topGo' GO categories
##
format(go_mf_topgenes[,c('associated_gene_name', 'logFC', 'FDR', 'fpkm_lps', 'fpkm_ctrl')], digits= 2)
dim(go_mf_topgenes)

barplot(sort(go_mf_topgenes$logFC), names.arg= go_mf_topgenes$associated_gene_name[order(go_mf_topgenes$logFC)], 
    las= 1, cex.names= 0.75, horiz= T, xlab= "Log2 fold change", cex.main=1,
    main= "Fold change for genes differentially expressed\nand mapped to GO molecular functions")
savePlot('M:/Documents/LabBook/LabBook_Figures/20100422_go_annotation_GOstats_Fig_4.emf', 'emf')

#barplot(t(as.matrix(log2(go_mf_topgenes[, c("fpkm_lps", "fpkm_ctrl")]))), 
#    beside= T, names.arg= go_mf_topgenes$associated_gene_name, cex.names= 0.75)

#heatmap(as.matrix(log2(go_mf_topgenes[, c("fpkm_lps", "fpkm_ctrl")])), 
#    labRow= go_mf_topgenes$associated_gene_name, labCol= c("LPS", "CTRL"), col= heat.colors(256), cexCol=1, scale=c("column"), Colv= NA, Rowv= NA, )

#______________________________________________________________________________
#
#                         BIOLOGICAL PROCESS
#______________________________________________________________________________

## Compile GOstats object
##
params <- GSEAGOHyperGParams(name = "S. scrofa custom GSEA based annot Params",
  geneSetCollection = gsc, geneIds = genes, universeGeneIds = universe,
  ontology = "BP", pvalueCutoff = 0.05, conditional = FALSE,
  testDirection = "over")

## Perform test(s)
## 
BP_over_condF <- hyperGTest(params)

# -------------------------------[ BP Visualization ]--------------------------

(test_summary<- summary(BP_over_condF, pvalue= 0.01, categorySize= 10))

htmlReport(BP_over_condF, file= 'bp_over_condF.html')


##
##  Plot GO terms
## 

top_n<- c(1:10) ## Lines in test_summary to plot in the graph
topGo<- test_summary$GOBPID[top_n] ## GO selected
gog<- GOGraph(topGo, GOBPPARENTS)  ## Prepare graph object
nodes_go<- nodes(gog) ## GO IDs in graph gog

## Convert GO accession to GO name
nodes_name<- vector(mode='character')
for(n in nodes_go){
    go_acc<- go_db[n]
    nodes_name<- append(nodes_name, Term(go_db[[n]]))
    }

## Make nodes redder if more significant (besed on pvalue)
alpha_gradient<- heat.colors(n=2, rep(seq(1, 0.3, length.out= 4), each= ceiling(length(top_n)/4)))[c(seq(1,10,by=2), seq(2,10, by=2))]
fillcol<- rep('lightblue', length(nodes_go)) ## Colour for GO terms not in the selected list (i.e. parents of thje selected terms)
fillcol[nodes_go %in% topGo]<- alpha_gradient

## see ?GraphvizAttributes for more attributes
node.attrs <-
    makeNodeAttrs(gog, shape="ellipse",
                  label = substr(nodes_name,1,12), ## Use substr() to shorten the names
                  fillcolor = fillcol,
                  fontsize= 24.0, ## Has no effect!
                  width= 2      ## Width of the nodes
                  ) 
windows(width= 17/2.54, height= 17/2.54)
plot(agopen(gog,
            recipEdges="distinct", layoutType="dot",
            nodeAttrs = node.attrs, name = "foo"))
# locator()
text(x= 0, y= 230, labels= "GO: Biological process\nTest for over-representation", adj=0)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100422_go_annotation_GOstats_Fig_2.jpeg', 'jpeg')


#______________________________________________________________________________
#
#                         CELLULAR COMPONENTS
#______________________________________________________________________________

## Compile GOstats object
##
params <- GSEAGOHyperGParams(name = "S. scrofa custom GSEA based annot Params",
  geneSetCollection = gsc, geneIds = genes, universeGeneIds = universe,
  ontology = "CC", pvalueCutoff = 0.05, conditional = FALSE,
  testDirection = "over")

## Perform test(s)
## 
CC_over_condF <- hyperGTest(params)

# -------------------------------[ CC Visualization ]--------------------------

(test_summary<- summary(CC_over_condF, pvalue= 0.01, categorySize= 10))

htmlReport(CC_over_condF, file= 'cc_over_condF.html')


##
##  Plot GO terms
## 

top_n<- c(1:10) ## Lines in test_summary to plot in the graph
topGo<- test_summary$GOCCID[top_n] ## GO selected
gog<- GOGraph(topGo, GOCCPARENTS)  ## Prepare graph object
nodes_go<- nodes(gog) ## GO IDs in graph gog

## Convert GO accession to GO name
nodes_name<- vector(mode='character')
for(n in nodes_go){
    go_acc<- go_db[n]
    nodes_name<- append(nodes_name, Term(go_db[[n]]))
    }

## Make nodes redder if more significant (besed on pvalue)
alpha_gradient<- heat.colors(n=2, rep(seq(1, 0.3, length.out= 4), each= ceiling(length(top_n)/4)))[c(seq(1,10,by=2), seq(2,10, by=2))]
fillcol<- rep('lightblue', length(nodes_go)) ## Colour for GO terms not in the selected list (i.e. parents of thje selected terms)
fillcol[nodes_go %in% topGo]<- alpha_gradient

## see ?GraphvizAttributes for more attributes
node.attrs <-
    makeNodeAttrs(gog, shape="ellipse",
                  label = substr(nodes_name,1,12), ## Use substr() to shorten the names
                  fillcolor = fillcol,
                  fontsize= 24.0, ## Has no effect!
                  width= 2      ## Width of the nodes
                  ) 
windows(width= 17/2.54, height= 17/2.54)
plot(agopen(gog,
            recipEdges="distinct", layoutType="dot",
            nodeAttrs = node.attrs, name = "foo"))
# locator()
text(x= 0, y= 1000, labels= "GO: Cellular component\nTest for over-representation", adj=0)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100422_go_annotation_GOstats_Fig_3.jpeg', 'jpeg')


# ---------------------------------[ Tritume ]---------------------------------

genes = c(1,10,100)
  evi =  c("TAS", "IPI", "IMP")## c("ND","IEA","IDA") ## c("IEA", "ISS", "IC") IDA ND  ## IEA ISS IC  TAS IPI NAS 
  GOIds = c("GO:0008150","GO:0008152","GO:0001666")
  frameData = data.frame(cbind(GOIds,evi,genes))

  library(AnnotationDbi)
  frame=GOFrame(frameData,organism="Homo sapiens")
  allFrame=GOAllFrame(frame)

  getGOFrameData(allFrame)
