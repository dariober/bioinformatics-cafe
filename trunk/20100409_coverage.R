# ---------------------------------------------------------------------
#  This script makes the barplots on LabBook 9 Apr 2010.
#
#  Input file has been prepared with 20100412_pileup_to_coverage.py
# ---------------------------------------------------------------------

cov_dt<- read.table('C:/Tritume/20100409_RNAseq_CTRL_sscrofa9.56.pileup.cov', sep='\t', header=TRUE)
names(cov_dt)

cov_cum<- tail(cumsum(rev(cov_dt$count_of_positions)), 10)

par(mfrow=c(1,2), mar=c(5,4,1,1), cex= 0.8, oma=c(0,0,3,0))
  barplot(cov_dt[1:10,2]/10^6, names.arg= cov_dt[1:10,1], xlab= 'Depth of coverage', 
    ylab= 'No. of positions (millions)', main= '')
  
  barplot(cov_cum/10^6, names.arg= 10:1, xlab= 'Cumulative coverage (from max to...)', ylab= '')

  legend('topleft', 
    legend= c(
      paste('n positions=', round(sum(cov_dt[,2])/10^6, 2), 'millions'),
      paste('max depth=', max(cov_dt[,1]), 'reads/base')
      ),
    bty = 'n', y.intersp= 1.5
    )
  title('Number of genome positions sequenced with depth x\n(library CTRL)', outer=T)
