## Description

Segment a sequence of observations in *states* using HMM model. Data can be modeled either as a discrete or as a normal distribution. 
Originally used to segment P-values along a chromsomes after having recoded p-values to discrete values. 
`segmentPvalueByHMM` is effectively a wrapper around the R package RHMM to handle data along chromosomes. 

## TODO

* Very important: When using normal distribution, the state ids are not unique as the counter is reset for every chunk of data 
processed.
