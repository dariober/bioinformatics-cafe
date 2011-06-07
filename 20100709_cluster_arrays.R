#
# This script performs some cluster analysis on the pig affymetrix arrays
# on BMDM.
# 

library(RODBC)
library(hopach)

# ---------------------------------[ Prepare array dataset ]------------------------------
#
# Read it if available:
# arr<- read.table('/exports/work/vet_roslin_nextgen/dario/R/20100709_hopach_microarray/array_bmdm_ct.txt', header= T, sep= '\t')

# Prepare it from CAGEdb:
#

conn_db<- odbcConnect(dsn= 'pgCAGE', case= 'nochange')

sqlQuery(conn_db,
  ## Generate crosstab table of all the arrays
  "SELECT cross_tab($$ 
                       SELECT 
                             affyarray_bmdm.probe_set_id,
                             CASE WHEN gene_symbol not like '---' 
                                  THEN gene_symbol || ' ' || affyarray_bmdm.probe_set_id
                                  ELSE affyarray_bmdm.probe_set_id END AS feature_id,
                             slide_id, 
                             intensity 
                       FROM affyarray_bmdm INNER JOIN porcine_na30_annot ON 
                            porcine_na30_annot.probe_set_id = affyarray_bmdm.probe_set_id
                     $$, 'tmp_array_ct'
                   );"
  )

arr<- sqlQuery(conn_db,
   ## Fetch crosstab table
   "SELECT * FROM tmp_array_ct ORDER BY probe_set_id;"
   )
arr[1:10,1:10]

## Write out crosstab table if required
write.table(arr, 'D:/Tritume/array_bmdm_ct.txt', sep= '\t', row.names= FALSE)

## Delete crosstab table from database and close connection
sqlQuery(conn_db, "drop table tmp_array_ct;")
odbcClose(conn_db)

# --------------------------[ Reformat array dataset ]--------------------------

## Remove columns of names (probe_set_id and feature_id) and assing rownames to feature_id
rownames(arr)<- arr[, "feature_id"]
probe_set_id<- as.vector(arr[, 'probe_set_id'])
arr<- arr[, 3:ncol(arr)]
arr[1:10,1:10]

## Log2 transformation
arr<- log2(arr)

# --------------------------[ hopach gene clustering ]--------------------------

## Calculate the variance of each probe
vars<- apply(arr, 1, var)
vars[1:10]

## Quantiles
quan<- quantile(vars)

## Upload an hopach object. Possibly created in Eddie and exported.
## load('U:/Documents/affimetrix/20100709_array_cluster_hopach/gene.hobj_1/gene.hobj.Rdata') ## This clustering using all probes with variance quantile >0.95 (1207 probes)
load('U:/Documents/affimetrix/20100709_array_cluster_hopach/gene.hobj_2/gene.hobj.Rdata')    ## This clustering using all probes with variance quantile >0.75 (6031 probes)


## Use probes in the top x quantile (e.g. 0.75 extracts the top 25% most variable probes)
## IMPORATANT: make sure the subset selected here agrees with the hopach object loaded above
q<- vars > quantile(vars, 0.75)
arr_quant<- arr[q, ]
probe_set_id<- probe_set_id[q]
dim(arr_quant)

## Calculate distance matrix (using distance metrix 'cosangle' as suggested in the vignettes)
gene.dist<- distancematrix(arr_quant, d= 'cosangle')
dim(gene.dist)
gene.dist[1:10,1:10]

##--------------- Skip this part if an hopach object is already loaded --------
##
## Function hopach() can take a while and it might be necessary to execute it on
## eddie, then save it with save()
## Run hopach

# gene.hobj<- hopach(arr_quant, dmat= gene.dist)

# -----------------------[ Cluster analysis ]----------------------------------

## Number of clusters
gene.hobj$clustering$k

gene.hobj$clustering$sizes

## It might be necessary to run this command on eddie. Then transfer files to a convenient dir (e.g. U:\Documents\affimetrix\20100709_array_cluster_hopach\hopach2tree)
hopach2tree(arr_quant, file= 'D:/Tritume/arr_quant.txt', 
  hopach.genes = gene.hobj, 
  # hopach.arrays = array.hobj, ## Uncomment this if array clustering has been done by hopach (see below)
  dist.genes = gene.dist)

## Link features to clusters
ann.clusters<- data.frame(cbind(
    feature_order= gene.hobj$final$order, 
    label= gene.hobj$cluster$labels), 
    label_final= gene.hobj$final$label,
    feature= rownames(arr_quant),
    probe_set_id= probe_set_id)

## Send clusters to CAGEdb
conn_db<- odbcConnect(dsn= 'pgCAGE', case= 'nochange')

sqlQuery(conn_db, 'DROP TABLE tmp_hopach')
sqlSave(conn_db, ann.clusters, tablename= 'tmp_hopach', rownames= F, colnames= F);
sqlQuery(conn_db, "COMMENT ON TABLE tmp_hopach IS 'Microarray probe clusters created in R by the hopach package. See script 20100709_cluster_arrays.R'")

ann.clusters[ann.clusters$feature == "Ssc.22002.2.A1_at", ]
ann.clusters[ann.clusters$feature == "TNF Ssc.100.1.S1_at", ] ## TNF

dim(ann.clusters[ann.clusters$label == 1200000,])

sort(unique(ann.clusters$label))
sort(unique(gene.hobj$final$label))[1:100]

# --------------------------[ Cluster ARRAYS ]----------------------------
#
# Using hclust in heatmap

## Correlation matrix
corr_mat<- cor(arr)

## Write out correlation matrix
write.table(file= 'D:/Tritume/array_corr.txt', x= corr_mat, sep= '\t')

## Heatmap of the arrays. This picture re-edited in Powerpoint.
heatmap(corr_mat, margins=c(10,10))


# --------------------------[ hopach array clustering ]------------------------

q<- vars > quantile(vars, 0.95)
arr_quant<- arr[q, ]
dim(arr_quant)

array.dist<- distancematrix(t(arr_quant), d= 'euclid')
array.hobj <- hopach(t(arr_quant), dmat= array.dist)
array.hobj$clust$k

# hopach2tree(arr_quant, file= 'D:/Tritume/arr_quant.arrays', hopach.genes = gene.hobj,
#   hopach.arrays = array.hobj, dist.arrays= array.dist)

par(mar= c(10,10,2,1))
dplot(array.dist, array.hobj, ord = "final", main = "", showclusters = TRUE, labels= colnames(arr_quant))

colnames(arr_quant)[array.hobj$final$order]

# ------------------------------[ Tritume ]------------------------------------

x<- 1:1000
y<- x + rnorm(n= 1000, mean= 10, sd= 100)
y[y<0]<- 0.1

cor(x,y)
cor(log2(x), log2(y))

hc.genedist<- dist()
hc<- hclust(dist(arr_quant))
plot(hc, xlim=c(0,1))
(rect.hclust(hc, k= 6))

km<- kmeans(arr_quant, 5)
plot(km)