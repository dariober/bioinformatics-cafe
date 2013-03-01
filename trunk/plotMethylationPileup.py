#!/usr/bin/env python

import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Plot the concatenated files from several outputs of mpileup2methylation.py
    
    Example input:
    --------------
    lambda_gi9626243        4       +       0       88      C       GSM882245.SRR449045
    lambda_gi9626243        7       +       1       87      C       GSM882245.SRR449045
    lambda_gi9626243        10      +       0       89      C       GSM882245.SRR449045
    lambda_gi9626243        11      +       0       89      C       GSM882245.SRR449045
    lambda_gi9626243        13      +       0       89      C       GSM882245.SRR449045

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--input', '-i',
                   required= True,
                   help='''Input file to plot.
                   ''')
parser.add_argument('--output', '-o',
                   required= True,
                   help='''Output file for the pdf file.
                   ''')
                   
args = parser.parse_args()

rscript= """
pile<- read.table('%(input)s', sep= '\t', stringsAsFactors= FALSE)
names(pile)<- c('chrom', 'pos', 'strand', 'cnt_met', 'cnt_unmet', 'base', 'library_id')
pile$depth<- pile$cnt_met + pile$cnt_unmet
pile$pct_met<- 100*(pile$cnt_met / pile$depth)

libraries<- sort(unique(pile$library_id))
pdf('%(output)s', height= 15/2.54, width= 22/2.54, pointsize= 9)
for(lib in libraries){
    pdata<- pile[which(pile$library_id == lib),]
    par(mfrow= c(2,1), las= 1, mgp= c(2, 0.5, 0), mar= c(1,4,3,1), oma= c(4, 1, 1, 2), xpd= NA)
    plot(pdata$pos / 1000, pdata$pct_met, ylim= c(0,100), xlim= c(1, max(pile$pos)/1000), type= 'h', col= 'grey', xlab= '', ylab= '%% methylated', bty= 'l', main= lib)
    abline(v= c(10, 20, 30, 38, 48), col= 'blue', lty= 'dotted', xpd= FALSE)
    arrows(x0= c(0, 20, 38), x1= c(10, 30, 48), y0= 102, y1= 102, lwd= 2, col= 'red', length= 0.05, angle= 90, code= 3)
    mtext(side= 3, line= 0.07, at= c(5, 25, 43), text= c('5mC', 'C', '5hmC'))

    plot(pdata$pos / 1000, pdata$depth+1, xlim= c(1, max(pile$pos)/1000), type= 'h', col= 'grey', xlab= 'Position (kb)', ylab= 'Read depth', bty= 'l', main= '', log= 'y')
    abline(v= c(10, 20, 30, 38, 48), col= 'blue', lty= 'dotted', xpd= FALSE)
    arrows(x0= c(0, 20, 38), x1= c(10, 30, 48), y0= max(pdata$depth), y1= max(pdata$depth), lwd= 2, col= 'red', length= 0.05, angle= 90, code= 3)
    mtext(side= 3, line= 0.07, at= c(5, 25, 43), text= c('5mC', 'C', '5hmC'))
    }
dev.off()
""" %{'input': args.input, 'output': args.output}

fout= open(args.input + '.R', 'w')
fout.write(rscript)
fout.close()

p= subprocess.Popen('Rscript %s ' %(args.input + '.R'), shell= True)
p.wait()
sys.exit()