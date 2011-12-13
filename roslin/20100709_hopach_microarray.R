# source("http://bioconductor.org/biocLite.R")
# install.packages(hopach)

## Run this script:
## source('/exports/work/vet_roslin_nextgen/dario/R/20100709_hopach_microarray/job_20100709_hopach_microarray.R')

library(hopach)

## Extract data from CAGE db with
# select cross_tab('select probe_set_id, slide_id, intensity from affyarray_bmdm', 'tmp_array_ct');
# copy tmp_array_ct to 'D:/Tritume/array_bmdm_ct.txt' with csv header delimiter E'\t';
#
# Move array_bmdm_ct.txt to /exports/work/vet_roslin_nextgen/dario/R/20100709_hopach_microarray/

arr<- read.table('/exports/work/vet_roslin_nextgen/dario/R/20100709_hopach_microarray/array_bmdm_ct.txt', header= T, sep= '\t')

## Remove first column (names) and assing rownames to probe IDs
rownames(arr)<- arr[,1]
arr<- arr[, 2:ncol(arr)]
arr[1:10,1:10]

## Log2 transformation
arr<- log2(arr)

vars<- apply(arr, 1, var)
vars[1:10]

## Quantiles
quan<- quantile(vars)

## Use probes in the top x quantile (e.g. 0.75 extracts the top 25% most variable probes)
q<- vars > quantile(vars, 0.75)
arr_quant<- arr[q, ]
dim(arr_quant)

## Calculate distance matrix (using distance metrix 'cosangle' as suggested in the vignettes)
gene.dist<- distancematrix(arr_quant, d= 'cosangle')
dim(gene.dist)
gene.dist[1:10,1:10]

## Run hopach (this takes a while)
gene.hobj<- hopach(arr_quant, dmat= gene.dist)

## Save gene.hobj. Use load() to reload it in a new session.
save(gene.hobj, file= '/exports/work/vet_roslin_nextgen/dario/R/20100709_hopach_microarray/gene.hobj.Rdata')

