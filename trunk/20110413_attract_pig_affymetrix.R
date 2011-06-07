library(attract)
library(porcine.db)
library(illuminaHumanv1.db)
library(RODBC)

# -----------------------------------------------------------------------------
# Some useful functions
# -----------------------------------------------------------------------------

synExpr2dataframe<- function(synset){
  ## Convert the slot "group" of a SynExpressionSet to dataframe with cluster
  ## slotNames(synset):
  ## [1] "groups"   "profiles"
  ##
  ## Use and ls() and get() the synset.
  
  myset.syn.dat<- as.data.frame(cbind(cluster_group= NA, probe_set_id= NA))[0,]
  for(i in 1:length(synset@groups)){
      cluster<- unlist(synset@groups[i])
      cluster_id<- rep(i, length(cluster))
      dat<- cbind(cluster_group= cluster_id, probe_set_id= cluster)
      myset.syn.dat<- rbind(myset.syn.dat, dat)
  }
  return(myset.syn.dat)
}

multi.synExpr2dataframe<- function(synsets){
  ## Convert an environment object with multiple synexpression groups to dataframe.
  ## In practice it applies synExpr2dataframe to each object:
  ## ls(synset)
  ## [1] "pway00480synexprs" "pway00590synexprs" "pway00830synexprs" ...
  synsets.df<- as.data.frame(cbind(keggid= NA, cluster_group= NA, probe_set_id= NA))[0,]
  all_keys<- ls(synsets)
  for(k in all_keys){
      syngroup<- get(k, synsets)
      keggid<- sub('pway', '', k); keggid<- sub('synexprs', '', keggid)
      kset<- synExpr2dataframe(syngroup)
      kset$keggid<- keggid
      kset<- kset[,c("keggid", "cluster_group", "probe_set_id")]
      synsets.df<- rbind(synsets.df, kset)
  }
  return(synsets.df)
}

pway2keggid<- function(pway){
    ## Convert pathway kegg id from attract to kegg id:
    ## "pway05416synexprs" >>> "05416"
    keggid<- sub('pway', '', k); keggid<- sub('synexprs', '', keggid)
}

row.aggr<- function(arow, aggr.by, fun){
    ## Applies a function to a row grouping by aggr.by
    x<- aggregate(arow, by= list(aggr.by), fun)$x
    return(x)
}

# -------------------------[ Import original dataset ]-------------------------

pig_array_ori<- read.table('U:/Documents/affymetrix/20110412_affymetrix_pig_rma.txt', sep='\t', header=TRUE, stringsAsFactors= FALSE)

colnames(pig_array_ori)[1:10]
rownames(pig_array_ori)[1:10]

## Input should look like this:
## --------------------------------------------------------------
# :> pig_array_ori[1:5,1:5]
#                s001_LWL1_0h_repl_1 s002_LWL1_2h_repl_1 s003_LWL1_7h_repl_1 s004_LWL1_24h_repl_1 s005_LWL1_0h_repl_2
# AFFX-BioB-3_at            8.249044            8.501810            8.496941             8.726919            8.400928
# AFFX-BioB-5_at            8.125517            8.366011            8.510154             8.691518            8.220821
# AFFX-BioB-M_at            8.381620            8.782618            8.844169             9.078234            8.695613
# AFFX-BioC-3_at            9.772719            9.957900            9.953191            10.206591            9.928164
# AFFX-BioC-5_at            9.087081            9.282562            9.338154             9.635377            9.253312
## --------------------------------------------------------------

## Average replicates
exprs.dat<- cbind(
    LWL1_0h=  rowMeans(pig_array_ori[, c(1, 5)]),      
    LWL1_2h=  rowMeans(pig_array_ori[, c(2, 6)]),
    LWL1_7h=  rowMeans(pig_array_ori[, c(3, 7)]),
    LWL1_24h= rowMeans(pig_array_ori[, c(4, 8)]),
    LWL4_0h=  rowMeans(pig_array_ori[, c(9, 13)]),
    LWL4_2h=  rowMeans(pig_array_ori[, c(10, 14)]),
    LWL4_7h=  rowMeans(pig_array_ori[, c(11, 15)]),
    LWL4_24h= rowMeans(pig_array_ori[, c(12, 16)]),
    LWL5_0h=  rowMeans(pig_array_ori[, c(17, 21)]),
    LWL5_2h=  rowMeans(pig_array_ori[, c(18, 22)]),
    LWL5_7h=  rowMeans(pig_array_ori[, c(19, 23)]),
    LWL5_24h= pig_array_ori[, c(20)]    ## Note: Array in column 25 "s024_LWL5_24h_repl_2" should be excluded (see LabBook 09/07/2010)
    )

## Re-order the array matrix to have time poimnts adjacent
exprs.dat<- exprs.dat[,c("LWL1_0h", "LWL4_0h", "LWL5_0h",
                         "LWL1_2h", "LWL4_2h", "LWL5_2h",
                         "LWL1_7h", "LWL4_7h", "LWL5_7h",
                         "LWL1_24h", "LWL4_24h", "LWL5_24h")
                      ]

exprs.dat[1:10,]

# -----------------------------------------------------------------------------
# Generate sample info
# -----------------------------------------------------------------------------

samp.info<- data.frame(array_id= colnames(exprs.dat),
  time_point= paste('tp_', sapply(strsplit(colnames(exprs.dat), '_'), '[', 2 ), sep= ''),
  pig_id= sapply(strsplit(colnames(exprs.dat), '_'), '[', 1 )
  )

## Sample data from the attract package
# data(exprs.dat)
# data(samp.info)

# -----------------------------------------------------------------------------
# Prepare expression set
# -----------------------------------------------------------------------------

loring.eset<- new("ExpressionSet")
loring.eset@assayData<- new.env()
assign('exprs', exprs.dat, loring.eset@assayData)
p.eset<- new("AnnotatedDataFrame", data= samp.info)
loring.eset@phenoData <- p.eset

# -----------------------------------------------------------------------------
# KEGG pathways showing the greatest discrimination between the different time
# points (or cell types tissue etc..) in the expression data set supplied.
# -----------------------------------------------------------------------------

## Put some large number for nperm ~ 1000
## See also gsealmPerm in GESAlm package
attractor.states<- findAttractors(loring.eset, "time_point", nperm= 1000,
    annotation= "porcine.db")
ranked_pathways<- attractor.states@rankedPathways
names(ranked_pathways)<- c('keggid', 'pathway', 'pvalue', 'size')
conn<- odbcConnect(dsn='pgVitelleschi')
sqlSave(conn, attractor.states@rankedPathways, tablename = 'attract_ranked_pathways', rownames= FALSE)
sqlQuery(conn, 'ALTER TABLE attract_ranked_pathways SET SCHEMA affymetrix')
sqlQuery(conn, "COMMENT ON TABLE attract_ranked_pathways IS 'Output of attract function findAttractors, slot rankedPathways. P-value for enrichment of kegg pathways. See 20110413_attract_pig_affymetrix.R' ")
odbcClose(conn)

## Test data
# attractor.states<- findAttractors(loring.eset, "celltype", nperm= 50,
#    annotation= "illuminaHumanv1.db")

slotNames(attractor.states)
## KEGG pathways ranked by pvalue
attractor.states@rankedPathways
  attractor.states@cellTypeTag

attractor.states@incidenceMatrix[, 1:10]

dim(attractor.states@incidenceMatrix)

# -----------------------------------------------------------------------------
# Remove uninformative genes
# -----------------------------------------------------------------------------

## Use output from LIMMA. This should come from 20110412_limma_affyarray_pig.R
remove.these.genes<- sqlQuery(conn, "
  select id, min(adjpval)
  from limma_toptable 
  group by id
  having min(adjpval) > 0.05;
")$id

length(remove.these.genes) ## 11614

#Consider uniformative probes which show no variation across the four time points
# remove.these.genes<- removeFlatGenes(loring.eset, "time_point", contrasts= NULL, limma.cutoff= 0.05)
# length(remove.these.genes)

# -----------------------------------------------------------------------------
# Find synexpression groups
# -----------------------------------------------------------------------------

# Identify transcriptionally coherent genes within a given pathway
# (synexpression groups).
# A synexpression group contains genes that share similar expression profiles
# across the tratment groups.
myset.syn<- findSynexprs(
    "00140", ## Pathway where to search for synexpression groups
    attractor.states, remove.these.genes)

## Get all groups where KEGG pvalue enrichments < p and Size large enough so that after
## removeing flat genes there is enough.
all_syn<- findSynexprs(
  subset(attractor.states@rankedPathways,  Pvalue < 0.05 & Size > 20)[,1],
  attractor.states, remove.these.genes)
slotNames(all_syn)

# -----------------------------------------------------------------------------
# INTERMEZZO: Get a specific KEGG group from the output of findSynexprs
# -----------------------------------------------------------------------------

## List groups in findSynexprs output:

# ls(all_syn) >>> "pway00480synexprs" "pway00590synexprs" "pway00830synexprs"  ...

## Get probes within clusters and expressions in a specific synexpression group:

#  get(ls(all_syn)[which(ls(all_syn) == "pway05320synexprs")], all_syn)
# An object of class "SynExpressionSet"
# Slot "groups":
#[[1]]
# [1] "Ssc.11025.1.S1_at"    "Ssc.11063.1.S1_at"    "Ssc.13780.1.S1_x_at"  "Ssc.13780.10.S1_x_at" "Ssc.13780.12.S1_a_at"
# [6] "Ssc.13780.12.S1_x_at" "Ssc.13780.3.S1_at"    "Ssc.13780.5.S1_x_at"  "Ssc.13780.6.S1_a_at"  "Ssc.13780.6.S1_x_at" 
#[11] "Ssc.13780.9.S1_a_at"  "Ssc.148.1.S1_at"      "Ssc.16160.1.S1_at"    "Ssc.222.1.S1_at"     
#
#[[2]]
#[1] "Ssc.14059.1.A1_at"   "Ssc.14497.1.S1_at"   "Ssc.15748.2.S2_at"   "Ssc.15748.3.S1_a_at" "Ssc.23660.1.S1_at"  
#[6] "Ssc.26880.1.A1_s_at" "Ssc.5915.1.S1_at"   
#
#
#Slot "profiles":
#       LWL1_0h   LWL4_0h  LWL5_0h  LWL1_2h   LWL4_2h  LWL5_2h   LWL1_7h  LWL4_7h   LWL5_7h  LWL1_24h  LWL4_24h  LWL5_24h
#[1,] 10.785343 10.371613 11.26561 11.21242 10.922053 11.62871 11.841710 11.73658 12.448339 12.141218 12.700364 12.753870
#[2,]  6.330959  6.341396  6.43962  7.94486  8.311066  8.75863  8.653192  8.79195  9.193986  7.102923  8.346249  7.798154
#
# -----------------------------------------------------------------------------
# Combine in a single dataframe the KEGG enrichment expression and the synexpression analysis.
# -----------------------------------------------------------------------------

names(attractor.states@rankedPathways)

all_synexpr<- data.frame( cbind(KEGGID= NULL, pway= NULL, pathway= NULL, pway_pvalue= NULL, group= NULL, exprs.dat[0,], probes= NULL ))
# i<- "pway00480synexprs"
# syn_profiles<- myset.syn@profiles
for( i in ls(all_syn) ){
    syn_profiles<- get( ls(all_syn)[which(ls(all_syn) == i)], all_syn )@profiles
    syn_groups<- get( ls(all_syn)[which(ls(all_syn) == i)], all_syn )@groups
    syn_groups<- sapply(syn_groups, paste, collapse= ', ')
    keggid<- rep(substring(i, 5, 9), nrow(syn_profiles) )
    cluster_group<- 1:nrow(syn_profiles)
    pathway<- as.vector(rep(attractor.states@rankedPathways[which(attractor.states@rankedPathways$KEGGID == keggid[1]), "PATHWAY"], nrow(syn_profiles)))
    pway_pvalue<- as.vector(rep(attractor.states@rankedPathways[which(attractor.states@rankedPathways$KEGGID == keggid[1]), "Pvalue"], nrow(syn_profiles)))
    all_synexpr<- rbind( all_synexpr, cbind(keggid, pway= rep(i, nrow(syn_profiles)), pathway= pathway, pway_pvalue= pway_pvalue, cluster_group, syn_profiles, syn_groups) )
}

write.table(all_synexpr, file= 'D:/Tritume/synexpression.txt', row.names= FALSE, col.names= TRUE, sep= '\t')

## Dataframe all_synexprs looks like this:
##
#:> all_synexpr[1:5,1:10]
#  keggid              pway                                      pathway       pway_pvalue cluster_group          LWL1_0h          LWL4_0h          LWL5_0h          LWL1_2h
#1  00480 pway00480synexprs                       Glutathione metabolism 0.004995004995005             1 10.7205609056249 10.5541796713826 10.7285279991508 10.6663232300038
#2  00480 pway00480synexprs                       Glutathione metabolism 0.004995004995005             2 5.92357105926505 5.86404040967929 5.90449619266517 5.93157681830835
#3  00590 pway00590synexprs                  Arachidonic acid metabolism 0.001998001998002             1 7.21613842396927 7.12989820782395  7.2827831637749 7.47991123099018
#4  00830 pway00830synexprs                           Retinol metabolism 0.005994005994006             1 4.66774441374573  4.4766880193078 4.97588756481018 4.76845107396061
#5  00980 pway00980synexprs Metabolism of xenobiotics by cytochrome P450 0.002997002997003             1 6.38529349306384   6.293435809267  6.5144083296087 6.41949787568575

sqlQuery(conn, "
    select read_table($$
    file: 'D:/Tritume/synexpression.txt',
    table: 'affymetrix.attract_synexpression',
    header: True,
    data_type: [['text']*3, 'double precision', 'int', ['double precision']*12, 'text'],
    overwrite: True $$)
")
sqlQuery(conn, "COMMENT ON TABLE affymetrix.attract_synexpression IS 'Output of R package attract, function findSynexprs. Probe clusters within signficantly enriched KEGG pathways. See 20110413_attract_pig_affymetrix.R' ")
## Now run 20110412_reformat_synexpression.sql to change the format of this table. It will produce the table attract_synexpression_probes
## Update comment for attract_synexpression_probes as appropriate
odbcClose(conn)

#-----------------------------[ Correlated probes ]-----------------------------

# dataframe with probes in kegg pathways
all_syn_df<- multi.synExpr2dataframe(all_syn)
all_syn_df[1:10,]

# group<- "pway00480synexprs"
all_syn_corrpartners<- as.data.frame(cbind(keggid= NA, cluster_group=NA, probe_set_id= NA))[0,]
for( group in ls(all_syn) ){
    syn_group<- get(group, all_syn)
    all_syn_corr<- findCorrPartners(syn_group, loring.eset, remove.these.genes, cor.cutoff = 0.9)
    all_syn_corr_df<- synExpr2dataframe(all_syn_corr)
    keggid<- sub('pway', '', group); keggid<- sub('synexprs', '', keggid);
    all_syn_corr_df$keggid<- keggid
    all_syn_corr_df<- all_syn_corr_df[, c("keggid", "cluster_group", "probe_set_id")]
    all_syn_corrpartners<- rbind(all_syn_corrpartners, all_syn_corr_df)
}
all_syn_corrpartners[1:10,]
dim(all_synexpr)

## Combine native probes in kegg pathways to the correlated one:
all_syn_corrpartners$kegg_probe<- 0
all_syn_df$kegg_probe<- 1 

sqlSave(conn, rbind(all_syn_df, all_syn_corrpartners), tablename= 'attract_synexpr_groups_kegg', rownames= FALSE )
sqlQuery(conn, 'ALTER TABLE attract_synexpr_groups_kegg SET SCHEMA affymetrix')
sqlQuery(conn, "COMMENT ON TABLE attract_synexpr_groups_kegg IS 'clusters produced by attract function findSynexprs() on the basis of kegg pathways. The output of findCorrPartners (probes correlated with r>0.9) concatenated. See 20110413_attract_pig_affymetrix.R'")
sqlQuery(conn, "COMMENT ON COLUMN attract_synexpr_groups_kegg.kegg_probe IS 'Does the probe belong to the current kegg pathway? If, 1 (Yes) the probe has been assigned by findSynexprs(), if 0 (No) the probe has been added by correlation using findCorrPartners(). See 20110413_attract_pig_affymetrix.R'")
# -----------------------------------------------------------------------------
# Plot clusters within the same pathway in the same graph
# -----------------------------------------------------------------------------

wdir<- 'F:/data/20110413_attract_pig_affymetrix'
dir.create(wdir)
wdir<- file.path(wdir, 'figure_misc')
dir.create(wdir)

setwd(wdir) 

## Pathways we want to plot:
pathways<- c("Antigen processing and presentation", "Cell adhesion molecules (CAMs)", "Chemokine signaling pathway",
             "Cytokine-cytokine receptor interaction", "Endocytosis", "Lysosome",
             "MAPK signaling pathway", "NOD-like receptor signaling pathway", "Phagosome",
             "RIG-I-like receptor signaling pathway", "TGF-beta signaling pathway", "Toll-like receptor signaling pathway")

## Abbreviations of pathway names for title of plots (order must match pathways vector!):
pathway_names<- c("Antigen proc. & pres.", "Cell adhesion mol.", "Chemokine sig. p'way",
             "Cytok.-cytok. rec. int.", "Endocytosis", "Lysosome",
             "MAPK sig. p'way", "NOD-like sig. p'way", "Phagosome",
             "RIG-I-like recept.", "TGF-beta sig. p'way", "TLR sig. p'way")

pathway_db<- as.data.frame(cbind(pathway= pathways, pathway_names), stringsAsFatctors= FALSE)

## Extract the selected pathways:
names(all_synexpr)
path<- all_synexpr[all_synexpr$pathway %in% pathways, c("pathway", "cluster_group", as.character(samp.info$array_id))]

## Average pigs and center on time point 0
ctrl_mean<- rowMeans( apply(path[, c("LWL1_0h", "LWL4_0h", "LWL5_0h")], 2, as.numeric) )  ## Generate a vector of the mean of the three 0 points for each cluster.
tp2mean_fc<- rowMeans( apply(path[, c("LWL1_2h", "LWL4_2h", "LWL5_2h")], 2, as.numeric) )
tp7mean_fc<- rowMeans( apply(path[, c("LWL1_7h", "LWL4_7h", "LWL5_7h")], 2, as.numeric) )
tp24mean_fc<- rowMeans( apply(path[, c("LWL1_24h", "LWL4_24h", "LWL5_24h")], 2, as.numeric) )

## Minimum and maximum at each time point
time_points<- rep(c(0,2,7,24), each= 3)
expr<- apply(path[, 3:ncol(path)], 2, as.numeric)
# path_max<- apply(expr, 1, row.aggr, time_points, function(x){mean(x) + sd(x)  } )
# path_min<- apply(expr, 1, row.aggr, time_points, function(x){mean(x) - sd(x)  } )
# path_max<- apply(expr, 1, row.aggr, time_points, function(x){mean(x) + (sd(x)/sqrt(length(x)))  } )
# path_min<- apply(expr, 1, row.aggr, time_points, function(x){mean(x) - (sd(x)/sqrt(length(x)))  } )
path_max<- apply(expr, 1, row.aggr, time_points, max )
path_min<- apply(expr, 1, row.aggr, time_points, min )
pway_fc<- data.frame(
                     pathway= path$pathway,
                     group= path$cluster_group,
                     tp_0h= 0,
                     tp_2h= tp2mean_fc - ctrl_mean,
                     tp_7h= tp7mean_fc - ctrl_mean,
                     tp_24h= tp24mean_fc - ctrl_mean,
                     tp_2h_top= path_max[2,] - ctrl_mean,
                     tp_7h_top= path_max[3,] - ctrl_mean,
                     tp_24h_top= path_max[4,] - ctrl_mean,
                     tp_2h_bottom= path_min[2,] - ctrl_mean,
                     tp_7h_bottom= path_min[3,] - ctrl_mean,
                     tp_24h_bottom= path_min[4,] - ctrl_mean
)

## Order by expression
pway_order<- data.frame(
    pathway= as.character(unique(pway_fc[order( -round(pway_fc$tp_2h, 1), -round(pway_fc$tp_7h, 1), -round(pway_fc$tp_24h, 1) ), "pathway"])),
    pway_order= 1:length(unique(pway_fc$pathway))
    )
pway_fc<- merge(pway_fc, pway_order, by= intersect("pathway", "pathway"))
pway_fc<- merge(pway_fc, pathway_db, by= intersect("pathway", "pathway"))
pway_fc<- pway_fc[order(pway_fc$pway_order), ]

## pathway<- pathways[1]

dev.off()
windows(width= 19/2.54, height= 13/2.54)
par(mfrow=c(3, 4), mar= c(0, 0.3, 1.5, 0), oma= c(1.5, 4, 1, 5), fg= "grey60", las= 1, cex= 0.65, bg= 'white')
# j<-1
pway_fc.dat<- apply(pway_fc[, 3:12], 2, as.numeric)
col<- c("firebrick4", "dodgerblue", "black", "darkolivegreen4")
pch<- c(19, 3, 21, 2)
for(j in 1:length(unique(pway_fc$pathway)) ){ 
  pathway<- as.character(unique(pway_fc$pathway)[j])
  pathway_name<- unique(pway_fc$pathway_name)[j]
  
  pway<- pway_fc[pway_fc$pathway == pathway, ]
  
  plot(x=c(0, 24), y= c(min(pway_fc.dat), max(pway_fc.dat)), type= 'n', main="", yaxt= 'n', xaxt= 'n', ylab= '', xlab= '', bty= 'o', col.axis= 'grey40')
  ## rect(xleft= par('usr')[1] + c(0, 6), ybottom= par('usr')[3], xright= par('usr')[1] + c(3, 9), ytop= par('usr')[4], col= 'grey95', border= 'transparent')
  rect(xleft= par('usr')[1], ybottom= par('usr')[4], xright= par('usr')[2], ytop= par('usr')[4] + strheight('text')*1.5, col= 'grey90', xpd= TRUE)
  box()
  for( i in 1:nrow(pway) ){
      polygon(x= c(0, 2, 7, 24, 24, 7, 2, 0),
              y= c(pway[i, c("tp_0h", "tp_2h_top", "tp_7h_top", "tp_24h_top")],
                   pway[i, c("tp_24h_bottom", "tp_7h_bottom", "tp_2h_bottom", "tp_0h")]),
              col= 'grey90', border= 'grey90'
      )
  }
  abline(h= 0, col= 'grey60', lwd= 1.25, lty= 'dotted')
    for( i in 1:nrow(pway) ){
      points(x=c(0, 2, 7, 24), y= pway[i, c("tp_0h", "tp_2h", "tp_7h", "tp_24h")], type= 'o', lwd= 1.75, main="", col= col[i], pch= pch[i], cex= 1)
    }
  if(j %in% c(9)){
        mtext(side= 1, text= c('0', '2', '7', '24 h'), at= c(0, 3, 7, 23), cex= 0.7, line= 0, col= 'grey25')        
  }
  if(j %in% c(1, 5, 9)){
        axis(side= 2)
  }
  mtext(side= 3, line= 0.2, pathway_name, adj= 0.5, col= 'grey25', cex= 0.7, font= 2)
  }
mtext(side= 2, text= 'log2 vs 0 h', cex= 0.9, line= 3, col= 'grey25', outer= TRUE, las= 0)
legend(x= 26.25, y= 7.5, legend= c('# 1', '# 2', '# 3', '# 4'), lwd= 1.2, col= col, pch= pch, xpd= NA, text.col= 'grey25')
savePlot(paste('all_clusters.avg.emf', sep= ''), 'emf')
savePlot('M:/Documents/LabBook/LabBook_Figures/20110413_all_clusters.avg.emf', 'emf')

# -----------------------------------------------------------------------------
# Cluster all the KEGG gene set
# -----------------------------------------------------------------------------

#
# Produce a dummy KEGG pathway containing all the probes and see what clusters
# are formed.
#

## Make a dummy AttractorModuleSet module:
attractor.states<- findAttractors(loring.eset, "time_point", nperm= 1,
    annotation= "porcine.db")

## Extend the incidenceMatrix to contain all the array probes.
## (However only those assigned to a KEGG pathway will be used).
new_mat<- matrix(data= 0, ncol= nrow(exprs.dat), nrow= nrow(attractor.states@incidenceMatrix),
    dimnames= list(rownames(attractor.states@incidenceMatrix), rownames(exprs.dat))
    )

## Produce custom gene set by assigning 1 to the probes to include:
new_mat["00140", ]<- 1

## Replace the dummy incidenceMatrix with the customized one:
attractor.states@incidenceMatrix<- new_mat

## Run findSynexprs
myset.syn<- findSynexprs(
    "00140", ## Pathway where to search for synexpression groups
    attractor.states, remove.these.genes)

# --------------------------[ Probes and Correlated probes ]-------------------

slotNames(myset.syn)
myset.syn@groups
corr_set_all<- findCorrPartners(myset.syn, loring.eset, remove.these.genes, cor.cutoff = 0.9)

myset.syn.dat<- synExpr2dataframe(corr_set_all)
myset.syn.dat$kegg_probe<- 0
myset.syn.df<- synExpr2dataframe(myset.syn)
myset.syn.df$kegg_probe<- 1

sqlSave(conn, rbind(myset.syn.df, myset.syn.dat), tablename= 'attract_synexpr_groups', rownames= FALSE )
sqlQuery(conn, 'ALTER TABLE attract_synexpr_groups SET SCHEMA affymetrix')
sqlQuery(conn, "COMMENT ON TABLE attract_synexpr_groups IS 'Clusters produced by attract function findSynexprs() from all the probes with kegg annotation grouped in the same dummy  kegg pathway. The output of findCorrPartners (probes correlated with r>0.9) concatenated. See 20110413_attract_pig_affymetrix.R'")
sqlQuery(conn, "COMMENT ON COLUMN attract_synexpr_groups.kegg_probe IS 'Does the probe belong to the current kegg pathway? If, 1 (Yes) the probe has been assigned by findSynexprs(), if 0 (No) the probe has been added by correlation using findCorrPartners(). See 20110413_attract_pig_affymetrix.R'")


head(myset.syn.dat)
head(myset.syn.df)

sapply(corr_set_all@groups, length)

## Convert log2 expressions to log2 FC relative to time-point 0
## I.e. substract from each pig-timepoint the time-point 0 for that pig.
pway_fc<- myset.syn@profiles
ctrl_tp<- pway_fc[, c('LWL1_0h', 'LWL4_0h', 'LWL5_0h')]
npigs<- length(unique(samp.info$pig_id))
for (i in seq(1, ncol(pway_fc), by= npigs) ){
    pway_fc[, i:(i+npigs-1)]<- pway_fc[, i:(i+npigs-1)] - ctrl_tp
    }
## Assign an ID to each profile, so it can be matched to the probe set after
## sorting.
rownames(pway_fc)<- 1:length(myset.syn@groups)

## Sort profile matrix by mean fold change:
## array_2h<- samp.info$array_id[samp.info$time_point == 'tp_2h']
## avg_2h<- rowMeans(pway_fc[, colnames(pway_fc) %in% array_2h])
avg_fc<- rowMeans(pway_fc)
pway_fc<- pway_fc[order(avg_fc, decreasing= TRUE), ]

## no. of probes in each cluster, sorted:
nprobes_sorted<- sapply( myset.syn@groups, length)[order(avg_fc, decreasing= TRUE) ]

dev.off()
windows(width= 20/2.54, height= 15/2.54)
## Make mfrow large enough to hold all the clusters. Or plot less clusters (see nrow(myset.syn@profiles) in for loop)
## par(mfrow=c(3, 5), mar= c(0, 0.3, 1.5, 0), oma= c(1.5, 3, 2, 0.2), fg= "grey60", las= 1, cex= 0.65, bg= 'white')
par(mfrow=c(5, 7), mar= c(0, 0.3, 1.5, 0), oma= c(1.5, 3, 2, 0.2), fg= "grey60", las= 1, cex= 0.65, bg= 'white')
col<- c("firebrick4", "dodgerblue", "grey40")
pch<- c(19, 3, 21)
## i<- 1
## pig<- 'LWL1'
## Make mfrow large enough to hold all the clusters. Or plot less clusters (see nrow(myset.syn@profiles) in for loop)
for(i in 1:nrow(myset.syn@profiles) ){
## for(i in 1:15 ){
  plot(x=c(0, 24), y= c(min(pway_fc), max(pway_fc)), type= 'n', main="", yaxt= 'n', xaxt= 'n', ylab= '', xlab= '', bty= 'o', col.axis= 'grey40')
  rect(xleft= par('usr')[1], ybottom= par('usr')[4], xright= par('usr')[2], ytop= par('usr')[4] + strheight('text')*1.5, col= 'grey90', xpd= TRUE)
  box()
  abline(h= 0, col= 'grey60', lwd= 1.25, lty= 'dotted')
    for(j in 1:length(unique(samp.info$pig_id)) ){
         pig<- unique(samp.info$pig_id)[j]
         tp_array<- as.character(samp.info$array_id[samp.info$pig_id == pig])
         pig_profile<- pway_fc[i, colnames(pway_fc) %in% tp_array]
         points(x=c(0, 2, 7, 24), y= pig_profile, type= 'o', lwd= 1.5, main="", col= col[j], pch= pch[j], cex= 1)
    }
  if(i %in% c(29, 31, 33, 35)){
        mtext(side= 1, text= c('0', '2', '7', '24 h'), at= c(0, 3, 7, 23), cex= 0.7, line= 0, col= 'grey25')        
  }
  if(i %in% c(1, 8, 15, 22, 29)){
        axis(side= 2)
  }
  mtext(side= 3, line= 0.2, paste('Cluster #', rownames(pway_fc)[i], ' (', nprobes_sorted[i], ')', sep= '' ), adj= 0.5, col= 'grey25', cex= 0.7, font= 2)  
}
mtext(side= 2, text= 'log2 vs 0 h', cex= 0.9, line= 2, col= 'grey25', outer= TRUE, las= 0)
mtext(side= 3, text= 'Clusters from all the probes mapping to KEGG pathways', outer= TRUE, col= 'grey25', line= 0.5)
savePlot('M:/Documents/LabBook/LabBook_Figures/all_kegg_probes_clustered_n35.emf', 'emf') ## Use this name for the plot showing 15 clusters
savePlot('M:/Documents/LabBook/LabBook_Figures/all_kegg_probes_clustered_n15.emf', 'emf') ## Use this name for the plot showing 15 clusters


## Export to postgres.
## Convert AttractorModuleSet to table

kegg_cluster<- as.data.frame(cbind(
    cluster_id= rep(1:length(myset.syn@groups), times= sapply(myset.syn@groups, length) ),
    probe_set_id= unlist(myset.syn@groups)),
    stringsAsFactors= FALSE)

conn<- odbcConnect(dsn= 'pgVitelleschi')
sqlSave(conn, kegg_cluster, 'attract_all_kegg', rownames= FALSE)
sqlQuery(conn, 'alter table attract_all_kegg set schema affymetrix')
sqlQuery(conn, "comment on table attract_all_kegg is 'Clusters produced by attract from all the probes mapping to any KEGG pathway (i.e. as if all the patwhays were collapsed in one), see 20110303_attarct_pig_affymetrix.R'")
odbcClose(conn)

# -----------------------------------------------------------------------------
# Find probes correlated with synexpression groups
# -----------------------------------------------------------------------------

# Enrich the synexpression groups, which contain genes annotates with KEGG, 
# with probes highly correlated to them.

synCorr<- data.frame(cbind(pway= NULL, group= NULL, corr_probes= NULL), stringsAsFactors= FALSE)[0,]
for( pway in ls(all_syn) ){
    syngroup<- get(ls(all_syn)[which(ls(all_syn) == pway)], all_syn)
    corr_set<- findCorrPartners(syngroup, loring.eset, remove.these.genes)
    corr_probes<- sapply(corr_set@groups, paste, collapse= ', ')
    synCorr<- rbind(synCorr,
                    cbind(
                        pway= rep(pway, length(syngroup@groups)),
                        group= seq(1:length(syngroup@groups)),
                        corr_probes= corr_probes
                    )
              )
}

conn<- odbcConnect(dsn= 'pgVitelleschi')
sqlQuery(conn, "CREATE TABLE affymetrix.attract_syncorr (pway text, cluster_group int, corr_probes text)")
write.table(synCorr, file= 'D:/Tritume/synCorr.txt', row.names= FALSE, col.names= TRUE, sep= '\t')
sqlQuery(conn, "COPY affymetrix.attract_syncorr FROM 'D:/Tritume/synCorr.txt' WITH CSV HEADER DELIMITER E'\t'")
sqlQuery(conn, "COMMENT ON TABLE affymetrix.attract_syncorr IS 'Output of R package attract, function findCorrPartners. Probes correlated with those in the synexpression groups. See 20110303_attract_pig_affymetrix.R' ")
odbcClose(conn)

# -----------------------------------------------------------------------------
# GO terms associated with synepression gruops
# -----------------------------------------------------------------------------

mapk.func <- calcFuncSynexprs(mapk.syn, attractor.states, "BP",
    annotation= "porcine.db")


# -----------------------------------------------------------------------------
# TRITUME
# -----------------------------------------------------------------------------

## Make a dummy AttractorModuleSet module:
attractor.states<- findAttractors(loring.eset, "time_point", nperm= 1,
    annotation= "porcine.db")

## Extend the incidenceMatrix to contain all the array probes.
new_mat<- matrix(data= 0, ncol= nrow(exprs.dat), nrow= nrow(attractor.states@incidenceMatrix),
    dimnames= list(rownames(attractor.states@incidenceMatrix), rownames(exprs.dat))
    )

## Produce custom gene set by assigning 1 to the probes to include:
new_mat["00140", ]<- 1

## Replace the dummy incidenceMatrix with the customized one:
attractor.states@incidenceMatrix<- new_mat

## Run findSynexprs
myset.syn<- findSynexprs(
    "00140", ## Pathway where to search for synexpression groups
    attractor.states, remove.these.genes)

plot(as.vector(myset.syn@profiles))

par(mfrow=c(4,4)) 
pretty.col <- rainbow(15) 
for( i in 1:15 ){
	plotsynexprs(myset.syn, tickMarks=c(6, 28, 47, 60), tickLabels=c("0", "2", "7", "24"), vertLines=c(12.5, 43.5, 51.5), index=i, 
			main=paste("Synexpression Group ", i, sep=""), col=pretty.col[i])
 }

sum(sapply(myset.syn@groups, length))

new_mat[1:10, 95:110]
new_incMat[1:10,1:10]

rowSums(attractor.states@incidenceMatrix[1:10, ])

attractor.states@incidenceMatrix[1:10, 1:10]

attractor.states@incidenceMatrix<- new_mat
plot(1:4, pch= c(19, 3, 21, 2))

cluster.dat[1:10,]


tnf<- cluster.dat[cluster.dat$gene_symbol == 'TNF', match(as.vector(samp.info$array_id), colnames(cluster.dat))][1,]
tnf<- tnf[order(samp.info$array_id)]

plot(1:4, tnf, type= 'n', xaxt= 'n')
points(x= 1, tnf)
samp.info

colnames(tnf)

order(c(1.1, 1.2, 1.0))


## Import from  postgres if necessary.
all_synexpr<- sqlFetch(conn, 'attract_synexpression')
odbcClose(conn)

## One graph per synexpression group. Images will be damped in this dir
setwd('F:/data/20110303_attract_pig_affymetrix/cluster_profiles') 

## pathway<- "Apoptosis"
for(pathway in as.character(all_synexpr$pathway)){
  pway<- all_synexpr[all_synexpr$pathway == pathway, ]
  dev.off()
  windows(width= 20/2.54, height= 4/2.54)
  par(mfrow=c(1,4), mar= c(1,3,1,1), oma= c(1, 1, 1, 0), fg= "grey60", las= 1, cex= 0.65, bg= 'white')
  df<- matrix(apply(pway[, 6:17], 2, as.numeric), ncol= 12)
  ctrl_mean<- apply( matrix(df[,1:3], ncol= 3), 1, mean) ## Generate a vector of the mean of the three 0 points for each cluster.
  df<- df - rep(ctrl_mean, 4)  ## Substract the ctrl mean to each time point to have the log2 fold change relative to 0 hours
  for( i in 1:nrow(pway) ){
      plot(x=1:12, y= df[i, ], type= 'n', main="", xaxt= 'n', ylab= 'n', xlab= '', bty= 'o', col.axis= 'grey40')
      rect(xleft= par('usr')[1] + c(0, 6), ybottom= par('usr')[3], xright= par('usr')[1] + c(3, 9), ytop= par('usr')[4], col= 'grey95', border= 'transparent')
      box()
      abline(h= 0, col= 'grey60', lwd= 2, lty= 'solid')
      points(x=c(2, 5, 8, 11), y= df[i, seq(1,12, by= 3)], type= 'o', lwd= 1.5, main="", col= "firebrick4", pch= 19, cex= 1.5)
      points(x=c(2, 5, 8, 11), y= df[i, seq(2,12, by= 3)], type= 'o', lwd= 1.5, main="", col= "dodgerblue", pch= 3, cex= 1.5)
      points(x=c(2, 5, 8, 11), y= df[i, seq(3,12, by= 3)], type= 'o', lwd= 1.5, main="", col= "grey40", pch= 21, cex= 1.5)
      mtext(side= 1, text= c('0 h', '2 h', '7 h', '24 h'), at= c(2, 5, 8, 11), cex= 0.7, line= 0.25, col= 'grey25')
    }
    mtext(side=3, text= pathway, adj= 0.07, font= 2, cex= 0.75, line= -0.5, col= 'grey25', outer= TRUE)
    mtext(side= 2, text= 'log2 vs 0 h', cex= 0.7, line= 0, col= 'grey25', outer= TRUE, las= 0)
  savePlot(paste(pathway, '.emf', sep= ''), 'emf')
  }


# -----------------------------------------------------------------------------
# Cluster visualization
# -----------------------------------------------------------------------------
## Data imported to postgres as above.
## Parsed by 20110303_attract_pig_affymetrix.sql and now back here.
## Probes have been linked to genes to get the expression levels.
library(RODBC)
conn<- odbcConnect(dsn= 'pgVitelleschi')

## Expression data in cross-tab format
cluster.dat<- sqlQuery(conn, 'select * from tmp_cluster_ct', stringsAsFactors= FALSE)
odbcClose(conn)
cluster.dat[1:10,]

setwd('F:/data/20110303_attract_pig_affymetrix/cluster_dendrograms')

pathways<- c("Antigen processing and presentation", "Cell adhesion molecules (CAMs)", "Chemokine signaling pathway",
             "Cytokine-cytokine receptor interaction", "Endocytosis", "Lysosome",
             "MAPK signaling pathway", "NOD-like receptor signaling pathway", "Phagosome",
             "RIG-I-like receptor signaling pathway", "TGF-beta signaling pathway", "Toll-like receptor signaling pathway")
for (pathway in pathways){
    ## pathway<- 'Apoptosis'
    pathway_cluster<- cluster.dat[which(cluster.dat$pathway == pathway), ]
    groups<- unique(pathway_cluster$cluster_group)
    for (group in groups){
        ## group<- 2 
        cluster.mat<- apply(pathway_cluster[which(pathway_cluster$cluster_group == group) , 4:ncol(cluster.dat)], 2, as.numeric)
        row.names(cluster.mat)<- pathway_cluster[which(pathway_cluster$cluster_group == group) , "gene_symbol"]
        
        dist_mat<- dist(cluster.mat)
        gene_cluster<- hclust(dist_mat, method= 'complete')
        
        dev.off()
        windows(width= 6/2.54, height= 8/2.54)
        par(mar= c(0, 0.25, 0, 7), oma=c(0,0,3,0), cex= 0.75)
        plot(as.dendrogram(gene_cluster), horiz= TRUE, axes= FALSE)
        mtext(text= paste(pathway, '\ngroup ', group, sep= ''), side= 3, line= 0, outer= TRUE, cex= 0.90, font =2)
        savePlot(paste(pathway, '-', group, '.dendrogram.emf', sep= ''), 'emf')
    }
}

#
# As above, but for selected pathways & clusters, ordered by expression
#

## Pathways of interest:
pathways<- c("Antigen processing and presentation", "Cell adhesion molecules (CAMs)", "Chemokine signaling pathway",
             "Cytokine-cytokine receptor interaction", "Endocytosis", "Lysosome",
             "MAPK signaling pathway", "NOD-like receptor signaling pathway", "Phagosome",
             "RIG-I-like receptor signaling pathway", "TGF-beta signaling pathway", "Toll-like receptor signaling pathway")

path<- all_synexpr[all_synexpr$pathway %in% pathways, c("pathway", "group", as.character(samp.info$array_id))]

## Average pigs and and center on time point 0

ctrl_mean<- rowMeans( apply(path[, c("LWL1_0h", "LWL4_0h", "LWL5_0h")], 2, as.numeric) )  ## Generate a vector of the mean of the three 0 points for each cluster.
tp2mean_fc<- rowMeans( apply(path[, c("LWL1_2h", "LWL4_2h", "LWL5_2h")], 2, as.numeric) )
tp7mean_fc<- rowMeans( apply(path[, c("LWL1_7h", "LWL4_7h", "LWL5_7h")], 2, as.numeric) )
tp24mean_fc<- rowMeans( apply(path[, c("LWL1_24h", "LWL4_24h", "LWL5_24h")], 2, as.numeric) )

pway_fc<- data.frame(pathway= path$pathway, group= path$group, tp_0h= 0, tp_2h= tp2mean_fc - ctrl_mean, tp_7h= tp7mean_fc - ctrl_mean, tp_24h= tp24mean_fc - ctrl_mean)

## Order by expression
pway_fc<- pway_fc[order( -round(pway_fc$tp_2h, 1), -round(pway_fc$tp_7h, 1), -round(pway_fc$tp_24h, 1) ), ]

## and plot
dev.off()
windows(width= 20/2.54, height= 16/2.54)
par(mfrow=c(5, 5), mar= c(0, 0, 0, 0), oma= c(1.5, 4, 3, 1), fg= "grey60", las= 1, cex= 0.65, bg= 'white')

for( i in 1:nrow(pway_fc) ){
    if( i > 25 ){
        break
    }
    plot(x=1:12, y= 1:12, ylim= c(min(pway_fc[,3:6]), max(pway_fc[,3:6])), type= 'n', main="", yaxt= 'n', xaxt= 'n', ylab= '', xlab= '', bty= 'o', col.axis= 'grey40')
    rect(xleft= par('usr')[1] + c(0, 6), ybottom= par('usr')[3], xright= par('usr')[1] + c(3, 9), ytop= par('usr')[4], col= 'grey95', border= 'transparent')
    box()
    abline(h= 0, col= 'grey60', lwd= 2, lty= 'solid')
    points(x=c(2, 5, 8, 11), y= pway_fc[i, c("tp_0h", "tp_2h", "tp_7h", "tp_24h")], type= 'o', lwd= 1.5, main="", col= "firebrick4", pch= 19, cex= 1.5)
    if(i %in% 21:25){
        mtext(side= 1, text= c('0 h', '2 h', '7 h', '24 h'), at= c(2, 5, 8, 11), cex= 0.7, line= 0, col= 'grey25')        
    }
    if(i %in% seq(1, 25, by=5)){
        axis(side= 2)
    }
    ## Add plot title
    if( (mean(as.numeric(pway_fc[i, c("tp_0h", "tp_2h", "tp_7h", "tp_24h")])) < 0) == TRUE ){
       side<- 3
    } else {
       side<- 1
    }
    mtext(side= side, line= -1.5, paste(substring(pway_fc$pathway[i], 1, 20), '.', sep= ''), adj= 0.2, col= 'grey25', cex= 0.7, font= 2)
  }
mtext(side= 2, text= 'log2 vs 0 h', cex= 0.9, line= 3, col= 'grey25', outer= TRUE, las= 0)
mtext(side=3, text= "Selected pathways", adj= 0.5, font= 2, cex= 1, line= +0.5, col= 'grey25', outer= TRUE)
savePlot('all_clusters_separate.emf', 'emf')
