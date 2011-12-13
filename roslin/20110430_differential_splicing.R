library(rpileup)

## TNF
rname<- 7
from<- 27658518
to<- 27661277
## ACTB
rname<- 3
from<- 3868536
to<- 3871915
## 11:48023539-48342202
rname<- 11
from<- 48023539
to<- 48342202
## CD14
rname<- 2
from<- 129356197
to<- 129357848
## CD58
rname<- 4
from<- 108346795
to<- 108569574
## CIDEB
rname<- 7
from<- 81260966
to<- 81459574

piles<- import.pileup(
    file.bam= c('F:/data/bamfiles/20110202_rnaseq_am_ctrl/accepted_hits.bam',
                 'F:/data/bamfiles/20110202_rnaseq_am_lps/accepted_hits.bam',
                 'F:/data/bamfiles/20110202_rnaseq_bmdm_ctrl/accepted_hits.bam',
                 'F:/data/bamfiles/20110202_rnaseq_bmdm_lps/accepted_hits.bam'),
    rname= rname, from= from, to= to
)
dim(piles)
features<- get.ensembl.features('sus', rname= rname, from= from, to= to)

plot.pileup(piles, overplot= FALSE,
    pileup.names= c('am_ctrl', 'am_lps', 'bmdm_ctrl', 'bmdm_lps'),
    col= c('steelblue', 'salmon'), type= 'l', oma= c(4,15,20,1), scale= 'max')
plot.features(features, transcript.spacing= 0.5, nlines= 3, box.cex= 1.025)
plot.exons(features)

dim(piles)

