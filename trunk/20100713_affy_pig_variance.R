#
#  Detecting microarray probes showing large variation
#  across pigs between 0h and treated (7 hours).
#

# -------------------------[ Import original dataset ]-------------------------

pig_array_ori<- read.table('U:/Documents/affimetrix/lwl_bmdm_lps.txt', 
  sep='\t', header=TRUE, stringsAsFactors= FALSE)

## Remove gene name from probe_set_id to return to the original Affymetrix names
## It is assumed the probe name doesn't have the colon ':' character.
#
probe_set_id<- sapply(pig_array_ori[,1], function(x) strsplit(x, ':')[[1]][1])

## Note: First take log2 of the intensity and then average
pig_array_log<- data.frame(probe_set_id= probe_set_id, log2(pig_array_ori[, 2:ncol(pig_array_ori)]))

## Average replicates
pig_array<- cbind(
LWL1_0h=  rowMeans(pig_array_log[, c(2, 6)]),      
LWL1_2h=  rowMeans(pig_array_log[, c(3, 7)]),
LWL1_7h=  rowMeans(pig_array_log[, c(4, 8)]),
LWL1_24h= rowMeans(pig_array_log[, c(5, 9)]),
LWL4_0h=  rowMeans(pig_array_log[, c(10, 14)]),
LWL4_2h=  rowMeans(pig_array_log[, c(11, 15)]),
LWL4_7h=  rowMeans(pig_array_log[, c(12, 16)]),
LWL4_24h= rowMeans(pig_array_log[, c(13, 17)]),
LWL5_0h=  rowMeans(pig_array_log[, c(18, 22)]),
LWL5_2h=  rowMeans(pig_array_log[, c(19, 23)]),
LWL5_7h=  rowMeans(pig_array_log[, c(20, 24)]),
LWL5_24h= pig_array_log[, c(21)]    ## Note: Array in column 25 "s024_LWL5_24h_repl_2" should be excluded (see LabBook 09/07/2010)
    )
rownames(pig_array)<- probe_set_id
pig_array[1:10,]

# -------------------------[ Prepare design table ]---------------------------

design<- as.data.frame(cbind(
  ## For convenience, prepare a design dataframe
  ## IMPORATANT: order here the same as in the array matrix
  array_id= c('LWL1_0','LWL1_2','LWL1_7','LWL1_24','LWL4_0','LWL4_2','LWL4_7','LWL4_24','LWL5_0','LWL5_2','LWL5_7','LWL5_24'),
  pig_id= c('LWL1','LWL1','LWL1','LWL1','LWL4','LWL4','LWL4','LWL4','LWL5','LWL5','LWL5','LWL5'),
  timepoint= c('0h','2h','7h','24h','0h','2h','7h','24h','0h','2h','7h','24h')
  ))

## Rename columns in pig_array
colnames(pig_array)<- as.character(design$array_id)

# -------------------------[ Estimate differential variance ]------------------

tp0h<- which(design$timepoint == '0h')
tp7h<- which(design$timepoint == '7h')

## Calculate F-statistic 0h vs 7h

p.ftest<- vector()
ftest<- vector()
for(i in 1:nrow(pig_array)){
    p.ftest<- append(p.ftest,     var.test(pig_array[i, tp0h], pig_array[i, tp7h])$p.value)
    ftest<- append(ftest,         var.test(pig_array[i, tp0h], pig_array[i, tp7h])$statistic)
    }

## Coefficient of variation

cv<- function(x){
    return(sd(x)/mean(x)) 
    }

cv0h<- apply(pig_array[,tp0h], 1, cv)
cv7h<- apply(pig_array[,tp7h], 1, cv)
avg0h<- apply(pig_array[,tp0h], 1, mean)
avg7h<- apply(pig_array[,tp7h], 1, mean)
sd0h<- apply(pig_array[,tp0h], 1, sd)
sd7h<- apply(pig_array[,tp7h], 1, sd)

# ----------------------------------[ Plot results ]-------------------------

## Plot CV ctrl vs LPS 7h
plot(x= cv0h, y= cv7h, xlab= 'CTRL 0 hours', ylab= 'LPS 7 hours', main= 'Coeffient of variation\n(cv= sd/mean)')
points(x= cv0h[p.ftest<0.01], y= cv7h[p.ftest<0.01], pch= 16, col= 'red')
legend('topleft', legend= 'p.val < 0.01 (F-test)', pch= 16, col= 'red')
length(p.ftest[p.ftest<0.01])
abline(a= 0, b= 1, col= 'blue', lwd= 2)

savePlot('M:/Documents/LabBook/LabBook_Figures/20100713_affy_pig_variance_CV.jpeg', 'jpeg')

## Plot profiles for some variable probes

probes<- c('Ssc.10097.1.A1_at', 'Ssc.15927.1.S1_at', 'Ssc.15927.2.S1_at', 'Ssc.9467.1.S1_at')
par(mfrow= c(2,2), mar= c(2,2,2,2), oma= c(2,2,2,0))
for(probe_id in probes){
    plot(y= pig_array[probe_set_id == probe_id, which(design$pig_id == 'LWL1')], 
        x= c(0,2,7,24), type= 'b', lwd= 1, col= 'blue', xlab= '', ylab= '',
        ylim= c(0, max(pig_array[probe_set_id == probe_id, ])), xaxt= 'n')
    axis(side= 1, at= c(0,2,7,24), labels= c(0,2,7,24))
    legend('bottomleft',legend= probe_id)
    points(y= pig_array[probe_set_id == probe_id, which(design$pig_id == 'LWL4')], 
        x= c(0,2,7,24), type= 'b', lwd= 1, col= 'red')
    points(y= pig_array[probe_set_id == probe_id, which(design$pig_id == 'LWL5')], 
        x= c(0,2,7,24), type= 'b', lwd= 1, col= 'green')
    }
mtext(text= c('Time point', 'Time point'), at= c(0.25, 0.75), side= 1, outer= T)
mtext(text= c('Log2 intensity', 'Log2 intensity'), at= c(0.25, 0.75), side= 2, outer= T, line= 0.5)
legend('bottomright', legend= c('Pig LWL1', 'Pig LWL4', 'Pig LWL5'), pch= 16,
    col= c('blue', 'red', 'green'), bty= 'n')
mtext(outer= T, text= 'Examples of probes showing pig-dependent reaction to LPS', cex= 1.25, side= 3, line= -0.75)
savePlot('M:/Documents/LabBook/LabBook_Figures/20100713_affy_pig_variance_probes.emf', 'emf')

## ----------------------[ Write out results ]---------------------
pig_var<- cbind(probe_set_id, pig_array[, c(tp0h, tp7h)], p_ftest= p.ftest, ftest, cv0h, cv7h, 
    avg0h, avg7h, sd0h, sd7h, cv_ratio= cv7h/cv0h
pig_var[1:10,]

write.table(pig_var, row.names= F,
    'U:/Documents/affimetrix/20100713_affy_pig_variance/pig_variance0-7h.txt', sep='\t'
    )