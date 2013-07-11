# source('~/svn_checkout/bioinformatics-misc/glm.R')


library(foreach)
library(doMC)

trans.arcsine <- function(x){
    "Transforms percentage data to linear via arcsine
    "
    asin(sign(x) * sqrt(abs(x)))
}


glmForeach<- function(i, bsobj, cond.1, cond.2, bs, thres, glm.family= 'binomial', ...){
    "Function to be executed inside foreach loop and returning a vector of required
    items.
    i:
        Row index to extract data from bsobj. Typically coming from foreach iteration 
    cond.1, cond.2:
        Vector of column indexes to identofy libraries belonging to condition 1 and 2.
    bs:
            Vector of characters to classify each library as bisulfite of oxBS.
            length(unique(bs)) must be 2
    thres:
        Exclude loci where tot count is less then thres (return NA to estimate
        and pvalue)
    family:
            String for Family of probablility distribution. Default is binomial but consider
            quasibinomial or gaussian instead. Valid options: 'gaussian', ''binomia', 'quasibinomial';
            If 'gaussian', counts are converted to percentage and asin(sqrt) transformed.
    "
    cnt_T<- as.integer(bsobj@tot_reads[i,])
    cnt_M<- as.integer(bsobj@cnt_met[i,])
    cnt_m<- cnt_T - cnt_M

    cond1.cnt_M<- sum(cnt_M[cond.1])
    cond1.cnt_m<- sum(cnt_m[cond.1])
    cond2.cnt_M<- sum(cnt_M[cond.2])         
    cond2.cnt_m<- sum(cnt_m[cond.2])

    if(sum(cnt_T) < thres){
        e<- NA
        p<- NA
    } else if (glm.family %in% c('binomial', 'quasibinomial')){
        lrfit <- glm(cbind(cnt_M, cnt_m) ~ as.factor(bs), family = glm.family)
        e<- summary(lrfit)$coef[, "Estimate"][2]
        pv<- grep('Pr\\(>\\|t|z\\|\\)', colnames(summary(lrfit)$coef), perl= TRUE)
        p<- summary(lrfit)$coef[, pv][2]
    } else if(glm.family %in% c('gaussian')){
        y.pct<- cnt_M / cnt_T
        y.asqr<- trans.arcsine(y.pct)
        lrfit <- glm(y.asqr ~ as.factor(bs), family = glm.family)
        e<- summary(lrfit)$coef[, "Estimate"][2]
        pv<- grep('Pr\\(>\\|t|z\\|\\)', colnames(summary(lrfit)$coef), perl= TRUE)
        p<- summary(lrfit)$coef[, pv][2]               
    } else {
        stop(sprintf('Unsopported family function: "%s"', family))
    }

    x<- c(locus= row.names(bsobj@tot_reads)[i],
        cnt.cond1.M= cond1.cnt_M,
        cnt.cond1.m= cond1.cnt_m,
        cnt.cond2.M= cond2.cnt_M,
        cnt.cond2.m= cond2.cnt_m,
        estimate= as.numeric(e),
        pvalue= as.numeric(p))
    return(x)
}

glmBS.2<- function(bsobj, bs, contrast, thres= 10, family= 'binomial', nthreads= 1, adjust.method= NA, ...){
    "
    _____________________________ Description ______________________________
    Fit a glm with to the each experiment or locus.
    
    _______________________________ Arguments ______________________________
    bsobj:
            BSdata object from which cnt_met, tot_reads, library_id, locus, bs can
            be extrated.
    bs:
            Vector of characters to classify each library as bisulfite of oxBS.
            length(unique(bs)) must be 2
    contrast:
            The contrast to apply to the levels of bs. Must be of length 2
            and with the same strings as in bs. E.g. c('BS', 'oxBS') to compare BS-oxBS
            or c('redBS', 'BS') to compare redBS-BS
    family:
            See glmForeach
    thres:
            See glmForeach
    nthreads:
            Number of threads to run in parallel (int passed to registerDoMC())
    adjust.method:
            Arg passed to p.adjust
    ...:
        Further args to glm()
    _____________________________ Value ____________________________________
    Data frame with the following columns:
    locus:
            The locus tested
    [...]
    estimate, pvalue:
            Estimate effect of difference contrast[1] - contrast[2] and p-value for significance
            of difference from 0 as obtained from glm()
    fdr:
            pvalue adjusted by method 'fdr'.
    
    MEMO: R defualt is to set contrast in reveresed alphanumeric order.
    For levels (contrats arg) BS and oxBS, R will test oxBS-BS meaning that positive
    estimates have oxBS > BS
    
    ____________________________ See also __________________________________
    
    http://en.wikipedia.org/wiki/Generalized_linear_model
    http://www.google.co.uk/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0CEEQFjAC&url=http%3A%2F%2Fwww.stat.umn.edu%2Fgeyer%2F5931%2Fmle%2Fglm.pdf&ei=7d-xUc6TAvSl0wWdzoGYCQ&usg=AFQjCNHaB9LhsQZnKol2Fl7IwwEd5bongw&sig2=qaAmsLYBC2Gxux_3939a5Q&bvm=bv.47534661,d.d2k&cad=rja
    "
    if(length(contrast) != 2){
        stop(sprintf('Only two contrasts allowed. Found %s', length(contrast)))
    }
    if(! all(contrast %in% bs)){
        stop(sprintf("Contrast not in arg bs"))
    }
    registerDoMC(nthreads)

    glm.family<- eval(parse(text= family))

    N<- nrow(bsobj@tot_reads)

    cond.1<- which(bs == contrast[1])
    cond.2<- which(bs == contrast[2])
    glmOut<- foreach(i= 1:N, .combine= rbind) %dopar% {
        if(i %% 1000 == 0){
            cat(sprintf('Locus N: %s/%s\n', i, N))
        }
        x<- glmForeach(i= i, bsobj= bsobj, cond.1= cond.1, cond.2= cond.2, bs= bs, thres= thres, glm.family= family)
    }
    glmOut<- data.frame(glmOut, stringsAsFactors= FALSE)

    glmOut[,1]<- as.character(glmOut[,1])
    glmOut[,2]<- as.integer(glmOut[,2])
    glmOut[,3]<- as.integer(glmOut[,3])
    glmOut[,4]<- as.integer(glmOut[,4])
    glmOut[,5]<- as.integer(glmOut[,5])
    glmOut[,6]<- as.numeric(glmOut[,6])
    glmOut[,7]<- as.numeric(glmOut[,7])

    if (is.na(adjust.method) == TRUE){
        glmOut$fdr<- NA
    } else {
        glmOut$fdr<- p.adjust(glmOut$pvalue, method= adjust.method)
    }
    return(glmOut)
}

# -------------------------------------------------------------------------------
# DEPRECATED VERSIONS
# -------------------------------------------------------------------------------

#glmBS.2<- function(bsobj, bs, contrast, thres= 10, family= 'binomial', nthreads= 1, adjust.method= NA, ...){
#    "
#    _____________________________ Description ______________________________
#    Fit a glm with to the each experiment or locus.
#    
#    _______________________________ Arguments ______________________________
#    bsobj:
#            BSdata object from which cnt_met, tot_reads, library_id, locus, bs can
#            be extrated.
#    bs:
#            Vector of characters to classify each library as bisulfite of oxBS.
#            length(unique(bs)) must be 2
#    contrast:
#            The contrast to apply to the levels of bs. Must be of length 2
#            and with the same strings as in bs. E.g. c('BS', 'oxBS') to compare BS-oxBS
#            or c('redBS', 'BS') to compare redBS-BS
#    family:
#            String for Family of probablility distribution. Default is binomial but consider
#            quasibinomial or gaussian instead. Valid options: 'gaussian', ''binomia', 'quasibinomial';
#            If 'gaussian', counts are converted to percentage and asin(sqrt) transformed.
#    thres:
#            Don't analyze loci with less than these many reads in total. These loci
#            don't appear at all in the output (to be chnaged?)
#    nthreads:
#            Number of threads to run in parallel (int passed to registerDoMC())
#     adjust.method:
#            Arg passed to p.adjust
#    ...:
#        Further args to glm()
#    _____________________________ Value ____________________________________
#    Data frame with the following columns:
#    locus:
#            The locus tested
#    pct.cond1.M, pct.cond2.M:
#            Avg PERCENTAGE methylated reads in condition 1 and 2 respectively. 1 and 2 assign
#            according to `contrast`.
#            The average is calculated by summing the total reads and methylated reads
#            from the same condition and then dividing by the number of libs in the condition.
#            It is *not* the average of the individual percentages within a condition.
#    avg.cnt_M, avg.cnt_c:
#            Average count of methylated and unmethyated reads per lib, regardless of condition
#    estimate, pvalue:
#            Estimate effect of difference contrast[1] - contrast[2] and p-value for significance
#            of difference from 0 as obtained from glm()
#    fdr:
#            pvalue adjusted by method 'fdr'.
#    
#    MEMO: R defualt is to set contrast in reveresed alphanumeric order.
#    For levels (contrats arg) BS and oxBS, R will test oxBS-BS meaning that positive
#    estimates have oxBS > BS
#    
#    ____________________________ See also __________________________________
#    
#    http://en.wikipedia.org/wiki/Generalized_linear_model
#    http://www.google.co.uk/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0CEEQFjAC&url=http%3A%2F%2Fwww.stat.umn.edu%2Fgeyer%2F5931%2Fmle%2Fglm.pdf&ei=7d-xUc6TAvSl0wWdzoGYCQ&usg=AFQjCNHaB9LhsQZnKol2Fl7IwwEd5bongw&sig2=qaAmsLYBC2Gxux_3939a5Q&bvm=bv.47534661,d.d2k&cad=rja
#    "
#    if(length(contrast) != 2){
#        stop(sprintf('Only two contrasts allowed. Found %s', length(contrast)))
#    }
#    if(! all(contrast %in% bs)){
#        stop(sprintf("Contrast not in arg bs"))
#    }
#    registerDoMC(nthreads)
#
#    glm.family<- eval(parse(text= family))
#
#    N<- nrow(bsobj@tot_reads)
#
#    estimate<- as.numeric(rep(NA, N))
#    pvalue<- as.numeric(rep(NA, N))
#    cond1.cnt_M<- as.integer(rep(NA, N))
#    cond1.cnt_m<- as.integer(rep(NA, N))
#    cond2.cnt_M<- as.integer(rep(NA, N))
#    cond2.cnt_m<- as.integer(rep(NA, N))
#    #avg.cnt_M<- as.numeric(rep(NA, N))
#    #avg.cnt_c<- as.numeric(rep(NA, N))
#    #pct.cond1.M<- as.numeric(rep(NA, N))
#    #pct.cond2.M<- as.numeric(rep(NA, N))
#
#    cond.1<- which(bs == contrast[1])
#    cond.2<- which(bs == contrast[2])
#
##    steps<- seq(1, floor(N / 10000)) * 10000
##    steps<- c(steps, steps[length(steps)] + (N %% 10000))
##    glmOut<- foreach(i= 1:N, .combine= rbind) %dopar% { 
#
#    for(i in seq(1, N)){
#        cnt_T<- as.numeric(bsobj@tot_reads[i,])
#        cnt_M<- as.numeric(bsobj@cnt_met[i,])
#        cnt_m<- cnt_T - cnt_M
#
##        pct_cond1_M<- 100 * ( sum(cnt_M[cond.1]) / sum(cnt_T[cond.1]) )
##        pct_cond2_M<- 100 * ( sum(cnt_M[cond.2]) / sum(cnt_T[cond.2]) )
##        sum_cnt_M<- sum(cnt_M)
##        sum_cnt_c<- sum(cnt_c)
#        cond1.cnt_M[i]<- sum(cnt_M[cond.1])
#        cond1.cnt_m[i]<- sum(cnt_m[cond.1])
#        cond2.cnt_M[i]<- sum(cnt_M[cond.2])         
#        cond2.cnt_m[i]<- sum(cnt_m[cond.2])
##        pct.cond1.M[i]<- sum(cnt_M[cond.1]) / sum(cnt_T[cond.1]) 
##        pct.cond2.M[i]<- sum(cnt_M[cond.2]) / sum(cnt_T[cond.2]) 
##        avg.cnt_M[i]<- mean(cnt_M)
##        avg.cnt_c[i]<- mean(cnt_c)
#
#        if(sum(cnt_T) < thres){
#            e<- NA
#            p<- NA
#        } else if (family %in% c('binomial', 'quasibinomial')){
#            lrfit <- glm(cbind(cnt_M, cnt_m) ~ as.factor(bs), family = glm.family, ...)
#            e<- summary(lrfit)$coef[, "Estimate"][2]
#            pv<- grep('Pr\\(>\\|t|z\\|\\)', colnames(summary(lrfit)$coef), perl= TRUE)
#            p<- summary(lrfit)$coef[, pv][2]
#        } else if(family %in% c('gaussian')){
#            y.pct<- cnt_M / cnt_T
#            y.asqr<- trans.arcsine(y.pct)
#            lrfit <- glm(y.asqr ~ as.factor(bs), family = glm.family, ...)
#            e<- summary(lrfit)$coef[, "Estimate"][2]
#            pv<- grep('Pr\\(>\\|t|z\\|\\)', colnames(summary(lrfit)$coef), perl= TRUE)
#            p<- summary(lrfit)$coef[, pv][2]               
#        } else {
#            stop(sprintf('Unsopported family function: "%s"', family))
#        }
#        estimate[i]<- e
#        pvalue[i]<- p
#        if(i %% 1000 == 0){
#            cat(sprintf('Locus N: %s/%s\n', i, N))
#        }
##        x<- c(cond1.cnt_M, cond1.cnt_m, cond2.cnt_M, cond2.cnt_m, e, p)
#    }
##     glmOut<- as.data.frame(glmOut)
#    glmOut<- data.frame(locus= row.names(bsobj@tot_reads), cond1.cnt_M, cond1.cnt_m, cond2.cnt_M, cond2.cnt_m, estimate, pvalue)
##     names(glmOut)<- c('cnt.cond1.M', 'cnt.cond1.m', 'cnt.cond2.M', 'cnt.cond2.m', 'estimate', 'pvalue')
#    if (is.na(adjust.method) == TRUE){
#        glmOut$fdr<- NA
#    } else {
#        glmOut$fdr<- p.adjust(glmOut$pvalue, method= adjust.method)
#    }
##     glmOut$locus<- row.names(bsobj@tot_reads)
#    return(glmOut)
##    glmOut<- data.frame(locus= row.names(bsobj@tot_reads), pct.cond1.M, pct.cond2.M, avg.cnt_M, avg.cnt_c, estimate, pvalue)
##    glmOut$pct.diff<- glmOut$pct.cond1.M - glmOut$pct.cond2.M
#}
#
#glmBS<- function(bsobj= NULL, cnt_M= NULL, cnt_c= NULL, library_id= NULL, bs= NULL, locus= NULL, contrast, thres= 10, family= 'binomial', ...){
#"
#_____________________________ Description ______________________________
#Fit a glm with to the each experiment or locus.
#
#_______________________________ Arguments ______________________________
#bsobj:
#         BSdata object from which cnt_met, tot_reads, library_id, locus, bs can
#         be extrated.
#cnt_M:
#         Vector of count of Methylated reads
#cnt_c:
#         Vector of count of un-methylated reads
#library_id:
#         Vector of Identifiers of the libraries
#bs:
#         Vector of characters to classify each library as bisulfite of oxBS.
#         length(unique(bs)) must be 2
#locus:
#         Identifier of each locus to test
#contrast:
#         The contrast to apply to the levels of bs. Must be of length 2
#         and with the same strings as in bs. E.g. c('BS', 'oxBS') to compare BS-oxBS
#         or c('redBS', 'BS') to compare redBS-BS
#family:
#          String for Family of probablility distribution. Default is binomial but consider
#          quasibinomial or gaussian instead. Valid options: 'gaussian', ''binomia', 'quasibinomial';
#          If 'gaussian', counts are converted to percentage and asin(sqrt) transformed.
#thres:
#         Don't analyze loci with less than these many reads in total. These loci
#         don't appear at all in the output (to be chnaged?)
#
#_____________________________ Value ____________________________________
#Data frame with the following columns:
#locus:
#         The locus tested
#pct.cond1.M, pct.cond2.M:
#         Avg PERCENTAGE methylated reads in condition 1 and 2 respectively. 1 and 2 assign
#         according to `contrast`.
#         The average is calculated by summing the total reads and methylated reads
#         from the same condition and then dividing by the number of libs in the condition.
#         It is *not* the average of the individual percentages within a condition.
#avg.cnt_M, avg.cnt_c:
#         Average count of methylated and unmethyated reads per lib, regardless of condition
#estimate, pvalue:
#         Estimate effect of difference contrast[1] - contrast[2] and p-value for significance
#         of difference from 0 as obtained from glm()
#fdr:
#         pvalue adjusted by method 'fdr'.
#
#____________________________ See also __________________________________
#
#http://en.wikipedia.org/wiki/Generalized_linear_model
#http://www.google.co.uk/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0CEEQFjAC&url=http%3A%2F%2Fwww.stat.umn.edu%2Fgeyer%2F5931%2Fmle%2Fglm.pdf&ei=7d-xUc6TAvSl0wWdzoGYCQ&usg=AFQjCNHaB9LhsQZnKol2Fl7IwwEd5bongw&sig2=qaAmsLYBC2Gxux_3939a5Q&bvm=bv.47534661,d.d2k&cad=rja
#"
#    if(length(contrast) != 2){
#         stop(sprintf('Only two contrasts allowed. Found %s', length(contrast)))
#    }
#    if(is.null(bsobj)){
#          if(length(unique(c(length(cnt_M), length(cnt_c), length(library_id), length(bs), length(locus)))) != 1){
#               stop('Vectors: cnt_M, cnt_c, library_id, bs, locus must be of equal length')
#          }
#          bsframe<- data.frame( locus= as.factor(locus), cnt_M, cnt_c, library_id= as.factor(library_id), bs= bs)
#    } else if(class(bsobj) == 'BSdata'){
#          ## Convert BSdata obj to R data.frame. Shouldn't do this. Keep using ff data frame instead!
#          cnt_M<- c(as.matrix((as.data.frame(as.ffdf(bsobj@cnt_met)))))
#          cnt_T<- c(as.matrix((as.data.frame(as.ffdf(bsobj@tot_reads)))))     
#          bsframe<- data.frame(locus= rep(rownames(bsobj@cnt_met), ncol(bsobj@cnt_met)),
#               cnt_M=  cnt_M,        ## c(bsobj@cnt_met),
#               cnt_c= cnt_T - cnt_M, ## c(bsobj@tot_reads - bsobj@cnt_met),
#               library_id= rep(colnames(bsobj@cnt_met), each= nrow(bsobj@cnt_met))
#          )
#          bsframe<- merge(bsframe, bsobj@design[, c('library_id', 'bs')], by.x= c('library_id'), by.y= c('library_id'), sort= FALSE)
#    } else {
#         stop('Invalid input.')
#    } 
#    bsframe<- bsframe[which(bsframe$bs %in% contrast),]
#    bsframe$bs<- factor(bsframe$bs, levels= contrast)
#    if(length(unique(bsframe$bs)) != 2){
#        stop(sprintf('Invalid number of contrast levels found: %s', unique(bsframe$contrast)))
#    }
#    estimate<- as.numeric()
#    pvalue<- as.numeric()
#    avg.cnt_M<- as.numeric()
#    avg.cnt_c<- as.numeric()
#    pct.cond1.M<- as.numeric()
#    pct.cond2.M<- as.numeric()
#    n<- 1
#    for(i in unique(bsframe$locus)){
#         bsdf<- bsframe[which(bsframe$locus == i), ]
#         filter<- which(bsdf$bs == contrast[1])
#         pct<- ( sum(bsdf$cnt_M[filter]) / sum(bsdf$cnt_M[filter] + bsdf$cnt_c[filter]) ) * 100
#         pct.cond1.M<- append(pct.cond1.M, pct)
#         filter<- which(bsdf$bs == contrast[2])
#         pct<- ( sum(bsdf$cnt_M[filter]) / sum(bsdf$cnt_M[filter] + bsdf$cnt_c[filter]) ) * 100
#         pct.cond2.M<- append(pct.cond2.M, pct)
#         avg.cnt_M<- append(avg.cnt_M, mean(bsdf$cnt_M))
#         avg.cnt_c<- append(avg.cnt_c, mean(bsdf$cnt_c))
#         if(sum(bsdf$cnt_M, bsdf$cnt_c) < thres){
#             e<- NA
#             p<- NA
#         } else if(family %in% c('binomial', 'quasibinomial')){
#             lrfit <- glm(cbind(cnt_M, cnt_c) ~ bs, family = eval(parse(text= family)), data= bsdf, ...)
#             e<- summary(lrfit)$coef[, "Estimate"][2]
#             pv<- grep('Pr\\(>\\|t|z\\|\\)', colnames(summary(lrfit)$coef), perl= TRUE)
#             p<- summary(lrfit)$coef[, pv][2]
#         } else if(family %in% c('gaussian')){
#             y.pct<- bsdf$cnt_M / rowSums(bsdf[, c('cnt_M', 'cnt_c')])
#             y.asqr<- trans.arcsine(y.pct)
#             lrfit <- glm(y.asqr ~ bsdf$bs, family = eval(parse(text= family)), ...)
#             e<- summary(lrfit)$coef[, "Estimate"][2]
#             pv<- grep('Pr\\(>\\|t|z\\|\\)', colnames(summary(lrfit)$coef), perl= TRUE)
#             p<- summary(lrfit)$coef[, pv][2]               
#         } else {
#               stop(sprintf('Unsopported family function: "%s"', family))
#         }
#        estimate<- append(estimate, e)
#        pvalue<- append(pvalue, p)
#        n<- n + 1
#        if(n %% 250 == 0){
#            cat(sprintf('Locus N: %s\n', n))
#        }
#    }
#    glmOut<- data.frame(locus= unique(bsframe$locus), pct.cond1.M, pct.cond2.M, avg.cnt_M, avg.cnt_c, estimate, pvalue)
#    glmOut$fdr<- p.adjust(pvalue, method= 'fdr')
#    glmOut$pct.diff<- glmOut$pct.cond1.M - glmOut$pct.cond2.M
#    return(glmOut)
#}
#
#empfdr<- function(glmbs, q= 0.05){
#     "Empirical FDR for significance of difference btw conditions in glmBS 
#     glmbs:
#          Data frame returned by glmBS
#     q:
#          Quantile to extract from negative p-values
#     "
#     pneg<- glmbs$pvalue[which(glmbs$pct.diff < 0)] ## Pvalues of negative percentages
#     qneg<- quantile(pneg, 0.05)
#     ppos<- glmbs$pvalue[which(glmbs$pct.diff > 0)]
#     qpos<- quantile(ppos, 0.05)
#     ## Number of positive % above qneg
#     pqthres<- length(which(ppos < qneg))
#     ## Number of negative % above qpos
#     nqthres<- length(which(pneg < qpos))
#     ## Positions exceeding 
#}
