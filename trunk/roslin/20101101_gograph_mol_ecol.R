library(GOstats)
library(GSEABase)
library(Rgraphviz)
library(gplots)


go_db<- as.list(GOTERM)

go_terms<- read.table('F:/Data/20101101_gographs_mol_ecol/go_terms.csv', header=T, sep= ',')
go_bp<- as.vector(go_terms[go_terms$Ontology == 'BP', 'GO_ID'])
gog<- GOGraph(go, GOBPPARENTS)

nodes_go<- nodes(gog) ## GO IDs in graph gog

## Convert GO accession to GO name
nodes_name<- vector(mode='character')
for(n in nodes_go){
    go_acc<- go_db[n]
    nodes_name<- append(nodes_name, Term(go_db[[n]]))
    }

fillcol<- rep('lightblue', length(nodes_go)) ## Colour for GO terms not in the selected list (i.e. parents of thje selected terms)
fillcol[nodes_go %in% go_bp]<- 'red'

node.attrs <-
    makeNodeAttrs(gog, shape="ellipse",
                  label = '', ## substr(nodes_name,1,4), ## Use substr() to shorten the names
                  fillcolor = fillcol,
##                  fontsize= 32.0, ## Has no effect!
                  width= 10,      ## Width of the nodes
                  height= 10
                  )

plot(agopen(gog,
            recipEdges="distinct", layoutType="dot",
            ## nodeAttrs = node.attrs, 
            name = "foo"))
g1<- layoutGraph(gog)
g_info<- nodeRenderInfo(g1)
text(x= g_info$nodeX[1:3], y= g_info$nodeY[1:3], labels= c(1,2,3), cex= 5)

savePlot('F:/Data/20101101_gographs_mol_ecol/bp_graph.emf', 'emf')


graph.par(list(nodes = list(col = "darkgreen", lty = "dotted",
  lwd = 2, fontsize = 24)))
g1<- layoutGraph(gog)
renderGraph(g1)

##
top_n<- c(1:10) ## Lines in test_summary to plot in the graph
topGo<- test_summary$GOMFID[top_n] ## GO selected

library("GO.db")
 g1 <- oneGOGraph("GO:0003680", GOMFPARENTS)
 g2 <- oneGOGraph("GO:0003701", GOMFPARENTS)
 g3 <- join(g1, g2)
 g4 <- GOGraph(c("GO:0003680", "GO:0003701"), GOMFPARENTS)

if( require("Rgraphviz") && interactive() )
  plot(g4)
