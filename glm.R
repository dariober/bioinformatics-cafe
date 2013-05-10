glmBS<- function(bsobj= NULL, cnt_M= NULL, cnt_c= NULL, library_id= NULL, bs= NULL, locus= NULL, contrast, thres= 10, ...){
"
_____________________________ Description ______________________________
Fit a glm with binomial error family to the each experiment or locus.

_______________________________ Arguments ______________________________
bsobj:
         BSdata object from which cnt_met, tot_reads, library_id, locus, bs can
         be extrated.
cnt_M:
         Vector of count of Methylated reads
cnt_c:
         Vector of count of un-methylated reads
library_id:
         Vector of Identifiers of the libraries
bs:
         Vector of characters to classify each library as bisulfite of oxBS.
         length(unique(bs)) must be 2
locus:
         Identifier of each locus to test
contrast:
         The contrast to apply to the levels of bs. Must be of length 2
         and with the same strings as in bs. E.g. c('BS', 'oxBS') to compare BS-oxBS
         or c('redBS', 'BS') to compare redBS-BS
thres:
         Don't analyze loci with less than these many reads in total. These loci
         don't appear at all in the output (to be chnaged?)

_____________________________ Value ____________________________________
Data frame with the following columns:
locus:
         The locus tested
pct.cond1.M, pct.cond2.M:
         Avg percentage methylated reads in condition 1 and 2 respectively. 1 and 2 assign
         according to `contrast`.
         The average is calculated by summing the total reads and methylated reads
         from the same condition and then dividing by the number of libs in the condition.
         It is *not* the average of the individual percentages within a condition.
avg.cnt_M, avg.cnt_c:
         Average count of methylated and unmethyated reads per lib, regardless of condition
estimate, pvalue:
         Estimate effect of difference contrast[1] - contrast[2] and p-value for significance
         of difference from 0 as obtained from glm()
fdr:
         pvalue adjusted by method 'fdr'.
"
    if(length(contrast) != 2){
         stop(sprintf('Only two contrasts allowed. Found %s', length(contrast)))
    }

    if(is.null(bsobj)){
          if(length(unique(c(length(cnt_M), length(cnt_c), length(library_id), length(bs), length(locus)))) != 1){
               stop('Vectors: cnt_M, cnt_c, library_id, bs, locus must be of equal length')
          }
          bsframe<- data.frame( locus= as.factor(locus), cnt_M, cnt_c, library_id= as.factor(library_id), bs= bs)
    } else if(class(bsobj) == 'BSdata'){
          bsframe<- data.frame(locus= rep(rownames(bsobj@cnt_met), ncol(bsobj@cnt_met)),
               cnt_M= c(bsobj@cnt_met),
               cnt_c= c(bsobj@tot_reads - bsobj@cnt_met),
               library_id= rep(colnames(bsobj@cnt_met), each= nrow(bsobj@cnt_met))
          )
          bsframe<- merge(bsframe, bsdata@design[, c('library_id', 'bs')], by.x= c('library_id'), by.y= c('library_id'), sort= FALSE)
    } else {
         stop('Invalid input.')
    } 
    bsframe<- bsframe[which(bsframe$bs %in% contrast),]
    bsframe$bs<- factor(bsframe$bs, levels= contrast)
    if(length(unique(bsframe$bs)) != 2){
        stop(sprintf('Invalid number of contrast levels found: %s', unique(bsframe$contrast)))
    }
    estimate<- as.numeric()
    pvalue<- as.numeric()
    avg.cnt_M<- as.numeric()
    avg.cnt_c<- as.numeric()
    pct.cond1.M<- as.numeric()
    pct.cond2.M<- as.numeric()
    n<- 1
    for(i in unique(bsframe$locus)){
         bsdf<- bsframe[which(bsframe$locus == i), ]
         filter<- which(bsdf$bs == contrast[1])
         pct<- ( sum(bsdf$cnt_M[filter]) / sum(bsdf$cnt_M[filter] + bsdf$cnt_c[filter]) ) * 100
         pct.cond1.M<- append(pct.cond1.M, pct)
         filter<- which(bsdf$bs == contrast[2])
         pct<- ( sum(bsdf$cnt_M[filter]) / sum(bsdf$cnt_M[filter] + bsdf$cnt_c[filter]) ) * 100
         pct.cond2.M<- append(pct.cond2.M, pct)
         avg.cnt_M<- append(avg.cnt_M, mean(bsdf$cnt_M))
         avg.cnt_c<- append(avg.cnt_c, mean(bsdf$cnt_c))
         if(sum(bsdf$cnt_M, bsdf$cnt_c) < thres){
             e<- NA
             p<- NA
         } else{
             lrfit <- glm(cbind(cnt_M, cnt_c) ~ bs, family = binomial, data= bsdf, ...)
             e<- summary(lrfit)$coef[, "Estimate"][2]
             p<- summary(lrfit)$coef[, "Pr(>|z|)"][2]
         }
        estimate<- append(estimate, e)
        pvalue<- append(pvalue, p)
        n<- n + 1
        if(n %% 250 == 0){
            cat(sprintf('Locus N: %s\n', n))
        }
    }
    glmOut<- data.frame(locus= unique(bsframe$locus), pct.cond1.M, pct.cond2.M, avg.cnt_M, avg.cnt_c, estimate, pvalue)
    glmOut$fdr<- p.adjust(pvalue, method= 'fdr')
    glmOut$pct.diff<- glmOut$pct.cond1.M - glmOut$pct.cond2.M
    return(glmOut)
}
