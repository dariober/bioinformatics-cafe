## R Code for calculating Transcript per Million

tpm<- function(counts, gene_length){
    # Transcripts per Million. See also
    # https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
    # Comparing TPM across samples makes sense IF the total number of transcripts
    # stays the same. This condition may be far from real.
    # 
    stopifnot(length(counts) == length(gene_length))
    stopifnot(is.numeric(counts))
    stopifnot(gene_length > 0)
    rate<- counts / gene_length  # Normalize gene count by length
    denom<- sum(rate)            # Total transcription
    xtpm<- rate / denom * 1e6    # Proportion each gene contributes to tot transcription;
    return(xtpm)                 # times 1M to have it has "count of mRNAs if the total
                                 # n of transcripts in the cell was 1M."
}

TPM<- function(counts, gene_length){
    # Apply tpm to vector or to each column of a matrix of counts. 
    if(is.vector(counts)){
        xtpm<- tpm(counts, gene_length)
    } else {
        xtpm<- apply(counts, 2, tpm, gene_length= gene_length)
    }
    return(xtpm)
}

## TEST
## =============================================================================
glen<- c(1, 10, 1000, 10000)
cnt<- c(10, 100, 10000, 300000)

x<- TPM(cnt, glen)

stopifnot(sum(x) == 1e6)
stopifnot(length(unique(x[1:3])) == 1)
stopifnot(x[1]*3 == x[4])

x<- TPM(cbind(a= cnt, b= cnt, c= cnt*2), glen)
stopifnot(colSums(x) == 1e6)
