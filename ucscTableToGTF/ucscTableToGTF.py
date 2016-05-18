#!/usr/bin/env python

import argparse
import sys
import urllib2

VERSION= '0.1.0'

class ucscLine:
    def __init__(self, line):
        line= line.strip().split('\t')
        self.transcript_id= line[1]
        self.chrom= line[2]
        self.strand= line[3]
        self.txStart= int(line[4])
        self.txEnd= int(line[5])
        self.cdsStart= int(line[6])
        self.cdsEnd= int(line[7])
        self.exonCount= int(line[8])
        self.exonStarts= [int(x) for x in line[9].split(',') if x != '']
        self.exonEnds= [int(x) for x in line[10].split(',') if x != '']
        self.score= int(line[11])
        self.name2= line[12]
        self.exonFrames= [int(x) for x in line[15].split(',') if x != '']
        self.attrs= 'gene_id "' + self.transcript_id + '"; ' + 'transcript_id "' + self.transcript_id + '"; ' + 'gene_name "' + self.name2 + '"; '
    ## Remeber to be careful with 0 and 1-based indexes: 
    ## ucsc table has starts 0 based
    ## GTF has starts 1 based so = UCSC-start+1
    ## END positions are the the same in both
    
    def checkCodonLenThree(self, gtf_codon_list):
        xlen= 0
        for x in gtf_codon_list:
            xlen += ((x[4] - x[3]) + 1)
        if xlen != 3:
            sys.stderr.write('\n' + line + '\n')
            sys.stderr.write('Invalid codon length: ' + str(xlen) + '\n\n')
            sys.exit(1)

    def checkCorrectCoords(self, gtf_line_list):
        for x in gtf_line_list:
            if x[3] > x[4]:
                xout= [str(x) for x in gtf_line_list]
                sys.stderr.write('\n' + line + '\n')
                sys.stderr.write('\n' + '\n'.join(xout) + '\n')
                sys.stderr.write('Error: start > end: %s-%s\n' %(str(x[3]), str(x[4])))
                sys.exit(1)

    def hasCDS(self):
        for x in self.exonFrames:
            if x > -1:
                return True
        return False

    def getTSS(self):
        if self.strand == '+':
            x= [self.chrom, '.', 'TSS', self.txStart+1, self.txStart+1, '0.000000', self.strand, '.', self.attrs]
        elif self.strand == '-':
            x= [self.chrom, '.', 'TSS', self.txEnd, self.txEnd, '0.000000', self.strand, '.', self.attrs]
        else:
            x= []
        return [x]
    def getTES(self):
        if self.strand == '+':
            x= [self.chrom, '.', 'TES', self.txEnd, self.txEnd, '0.000000', self.strand, '.', self.attrs]
        elif self.strand == '-':
            x= [self.chrom, '.', 'TES', self.txStart+1, self.txStart+1, '0.000000', self.strand, '.', self.attrs]
        else:
            x= []
        return [x]
    def getGtfExons(self):
        exons= []
        for i in range(len(self.exonStarts)):
            exon= [self.chrom, '.', 'exon', self.exonStarts[i]+1, self.exonEnds[i], '0.000000', self.strand, '.', self.attrs]
            exons.append(exon)
        self.checkCorrectCoords(exons)
        return exons

    def getGtfIntrons(self):
        introns= []
        for i in range(self.exonCount-1):
            introns.append( [self.chrom, '.', 'intron', self.exonEnds[i]+1, self.exonStarts[i+1], '0.000000', self.strand, '.', self.attrs] )
        self.checkCorrectCoords(introns)
        return introns 

    def get_stop_codon(self):
        stop_codons= [] ## You have only one stop per transcript but it can be split across to exons
        if len(set(self.exonFrames)) == 1 and self.exonFrames[0] == -1:
            return []
        if self.strand == '+':
            codonEnd= self.cdsEnd     ## This is correct
            codonStart= self.cdsEnd-2 ## This has to be checked
            ## You have to account for the case where the extended codon spills over the previous exon:
            """
                   Z    ZZ
              EEEEEE    EEEEEEEEE
            """
            # codon_exon -> the exon containing the codonEnd
            codon_exon_idx= None
            for i in range(self.exonCount):
                s, e= self.exonStarts[i]+1, self.exonEnds[i]
                if s <= codonEnd <= e:
                    codon_exon_idx= i
                    break
            
            if codon_exon_idx is not None:
                spill= (codonEnd - 2) - (self.exonStarts[codon_exon_idx] + 1)
                if spill < 0:
                    # if spill < 0, the codon spans two exons. Get end of previous exon 
                    # and set codonStart to fill up the spill
                    codonStart= (self.exonEnds[codon_exon_idx - 1] + 1) + spill ## NB spill is -ve so you are moving left with +spill!
                    ## Now you need to split the codon in two chunks: 
                    ## codonStart->previousExonEnd
                    stop_codons.append([self.chrom, '.', 'stop_codon', codonStart, self.exonEnds[codon_exon_idx - 1], '0.000000', self.strand, '.', self.attrs])
                    ## codingExonStart->codonEnd
                    stop_codons.append([self.chrom, '.', 'stop_codon', self.exonStarts[codon_exon_idx]+1, codonEnd, '0.000000', self.strand, '.', self.attrs])
                else:
                    stop_codons.append([self.chrom, '.', 'stop_codon', codonStart, codonEnd, '0.000000', self.strand, '.', self.attrs])    
            else:
                stop_codons.append([self.chrom, '.', 'stop_codon', codonStart, codonEnd, '0.000000', self.strand, '.', self.attrs])

        elif self.strand == '-':
            codonStart= self.cdsStart+1 ## This is correct
            codonEnd= self.cdsStart+3
            ## You have to account for the case where the extended codon spills over the previous exon:
            """
                   Z    ZZ
              EEEEEE    EEEEEEEEE
            """
            # codon_exon -> the exon containing the codonEnd
            codon_exon_idx= None
            for i in range(self.exonCount):
                s, e= self.exonStarts[i]+1, self.exonEnds[i]
                if s <= codonStart <= e:
                    codon_exon_idx= i
                    break
            if codon_exon_idx is not None:
                spill= (codonStart+2) - self.exonEnds[codon_exon_idx]
                if spill > 0:
                    ## First chunk
                    stop_codons.append([self.chrom, '.', 'stop_codon', codonStart, self.exonEnds[codon_exon_idx], '0.000000', self.strand, '.', self.attrs])
                    ## Second chunk
                    codonEnd= self.exonStarts[codon_exon_idx + 1] + spill
                    stop_codons.append([self.chrom, '.', 'stop_codon', self.exonStarts[codon_exon_idx+1]+1, codonEnd, '0.000000', self.strand, '.', self.attrs])
                else:
                    stop_codons.append([self.chrom, '.', 'stop_codon', codonStart, codonEnd, '0.000000', self.strand, '.', self.attrs])
            else:
                stop_codons.append([self.chrom, '.', 'stop_codon', codonStart, codonEnd, '0.000000', self.strand, '.', self.attrs])            
        else:
            x= []
        self.checkCorrectCoords(stop_codons)
        self.checkCodonLenThree(stop_codons)
        return stop_codons

    def get_start_codon(self):
        start_codons= [] ## You have only one stop per transcript but it can be split across to exons
        if len(set(self.exonFrames)) == 1 and self.exonFrames[0] == -1:
            return []
        if self.strand == '+':
            codonStart= self.cdsStart+1     ## This is correct
            codonEnd= self.cdsStart+3 ## This has to be checked
            ## You have to account for the case where the extended codon spills over the previous exon:
            """
                   A    AA
              EEEEEE    EEEEEEEEE
            """
            # codon_exon -> the exon containing the codonEnd
            codon_exon_idx= None
            for i in range(self.exonCount):
                s, e= self.exonStarts[i]+1, self.exonEnds[i]
                if s <= codonStart <= e:
                    codon_exon_idx= i
                    break
            if codon_exon_idx is not None:
                spill= (codonStart+2) - self.exonEnds[codon_exon_idx]
                if spill > 0:
                    ## First chunk
                    start_codons.append([self.chrom, '.', 'start_codon', codonStart, self.exonEnds[codon_exon_idx], '0.000000', self.strand, '.', self.attrs])
                    ## Second chunk
                    codonEnd= self.exonStarts[codon_exon_idx + 1] + spill
                    start_codons.append([self.chrom, '.', 'start_codon', self.exonStarts[codon_exon_idx+1]+1, codonEnd, '0.000000', self.strand, '.', self.attrs])
                else:
                    start_codons.append([self.chrom, '.', 'start_codon', codonStart, codonEnd, '0.000000', self.strand, '.', self.attrs])
            else:
                start_codons.append([self.chrom, '.', 'start_codon', codonStart, codonEnd, '0.000000', self.strand, '.', self.attrs])

        if self.strand == '-':
            codonEnd= self.cdsEnd     ## This is correct
            codonStart= self.cdsEnd-2 ## This has to be checked
            """
                   A    AA
              EEEEEE    EEEEEEEEE
            """
            # codon_exon -> the exon containing the codonEnd
            codon_exon_idx= None
            for i in range(self.exonCount):
                s, e= self.exonStarts[i]+1, self.exonEnds[i]
                if s <= codonEnd <= e:
                    codon_exon_idx= i
                    break
            
            if codon_exon_idx is not None:
                spill= (codonEnd - 2) - (self.exonStarts[codon_exon_idx] + 1)
                if spill < 0:
                    # if spill < 0, the codon spans two exons. Get end of previous exon 
                    # and set codonStart to fill up the spill
                    codonStart= (self.exonEnds[codon_exon_idx - 1] + 1) + spill ## NB spill is -ve so you are moving left with +spill!
                    ## Now you need to split the codon in two chunks: 
                    ## codonStart->previousExonEnd
                    start_codons.append([self.chrom, '.', 'start_codon', codonStart, self.exonEnds[codon_exon_idx - 1], '0.000000', self.strand, '.', self.attrs])
                    ## codingExonStart->codonEnd
                    start_codons.append([self.chrom, '.', 'start_codon', self.exonStarts[codon_exon_idx]+1, codonEnd, '0.000000', self.strand, '.', self.attrs])
                else:
                    start_codons.append([self.chrom, '.', 'start_codon', codonStart, codonEnd, '0.000000', self.strand, '.', self.attrs])    
            else:
                start_codons.append([self.chrom, '.', 'start_codon', codonStart, codonEnd, '0.000000', self.strand, '.', self.attrs])
        self.checkCorrectCoords(start_codons)
        self.checkCodonLenThree(start_codons)
        return start_codons

    def getGtfCDS(self):
        if len(set(self.exonFrames)) == 1 and self.exonFrames[0] == -1:
            return []
        CDSs= []
        for i in range(len(self.exonStarts)):
            exStart= self.exonStarts[i]
            exEnd= self.exonEnds[i]
            if self.exonFrames[i] == 0:
                frame= 0
            elif self.exonFrames[i] == 1:
                frame= 2
            elif self.exonFrames[i] == 2:
                frame= 1
            elif self.exonFrames[i] == -1:
                continue
            else:
                sys.stderr.write('\nUnexpected frame: %s \n' %(frame))
                sys.exit(1)

            if self.strand == '+':
                if (exStart <= self.cdsStart and exEnd >= self.cdsStart):
                    ## >>>>cdsStart>>>>
                    CDSs.append( [self.chrom, '.', 'CDS', self.cdsStart + 1, exEnd, '0.000000', self.strand, frame, self.attrs] )
                elif (exStart >= self.cdsStart and exEnd <= self.cdsEnd): 
                    ## cdsStart >>>>>>>> cdsEnd 
                    CDSs.append( [self.chrom, '.', 'CDS', exStart  + 1, exEnd, '0.000000', self.strand, frame, self.attrs] )
                elif (exStart <= self.cdsEnd and exEnd >= self.cdsEnd):
                    ## >>>>>>>>cdsEnd>>>>>>
                    CDSs.append( [self.chrom, '.', 'CDS', exStart + 1, self.cdsEnd - 3, '0.000000', self.strand, frame, self.attrs] )
                else:
                    continue
            elif self.strand == '-':
                if (exStart <= self.cdsStart and exEnd >= self.cdsStart):
                    CDSs.append( [self.chrom, '.', 'CDS', self.cdsStart + 1 + 3, exEnd, '0.000000', self.strand, frame, self.attrs] )
                elif (exStart >= self.cdsStart and exEnd <= self.cdsEnd): 
                    CDSs.append( [self.chrom, '.', 'CDS', exStart  + 1, exEnd, '0.000000', self.strand, frame, self.attrs] )
                elif (exStart <= self.cdsEnd and exEnd >= self.cdsEnd):
                    CDSs.append( [self.chrom, '.', 'CDS', exStart + 1, self.cdsEnd, '0.000000', self.strand, frame, self.attrs] )
                else:
                    continue
            else:
                continue
        ## Remove CDS where start > end. This can happen if the start/stop codon is split between exons
        CDSs= [x for x in CDSs if x[3] <= x[4]]
        
        ## Correct the end position of the last CDS by clipping the stop_codon away and anything afterwards
        if self.strand == '+':
            """
            CCCCCCCC...
                    ZZZ
            """
            start_of_stop_cod= self.get_stop_codon()[0][3]
            end_of_last_cds= CDSs[-1][4]
            if end_of_last_cds >= start_of_stop_cod:
                CDSs[-1][4]= start_of_stop_cod - 1
            """
            . is the bit to clip
                        ZZZ <- Stop codon
            CCCCCCCCCCCC... <- Last CDS
            """
        if self.strand == '-':
            """
                            aaa     <- Start codon
            ccccccccccccccccccc.... <- Rightmost CDS
                              *     <- end_of_start_cod
                                  * <- end_of_last_cds
            """
            end_of_start_cod= self.get_start_codon()[-1][4]
            end_of_last_cds= CDSs[-1][4]
            if end_of_last_cds > end_of_start_cod:
                CDSs[-1][4]= end_of_start_cod 
            """
            zzz
            ...cccccccc <- First CDS
            """
            end_of_stop_cod= self.get_stop_codon()[-1][4]
            start_of_first_cds= CDSs[0][3]
            if start_of_first_cds <= end_of_stop_cod:
                CDSs[0][3]= end_of_stop_cod + 1
        CDSs= [x for x in CDSs if x[3] <= x[4]] ## Again, remove possibly wring assignments
        self.checkCorrectCoords(CDSs)
        return CDSs
    

    def getGtf_utr_forward(self):
        utrs= []
        if self.strand == '-':
            return utrs
        if not self.hasCDS():
            for s, e in zip(self.exonStarts, self.exonEnds):
                xutr= utrs.append( [self.chrom, '.', '5utr', s+1, e, '0.000000', self.strand, -1, self.attrs] )
            self.checkCorrectCoords(utrs)
            return utrs            
        
        end_of_5utr= self.get_start_codon()[0][3] - 1
        start_of_3utr= self.get_stop_codon()[-1][4] + 1

        for i in range(self.exonCount):
            s, e= self.exonStarts[i]+1, self.exonEnds[i]
            if e <= end_of_5utr:
                ## Exon is fully 5utr
                xutr= utrs.append( [self.chrom, '.', '5utr', s, e, '0.000000', self.strand, -1, self.attrs] )
            if s <= end_of_5utr and e >= end_of_5utr:
                ## Exon is partially a 5utr
                xutr= utrs.append( [self.chrom, '.', '5utr', s, end_of_5utr, '0.000000', self.strand, -1, self.attrs] )
            if s <= start_of_3utr and e >= start_of_3utr:
                ## Exon is partially a 3utr
                xutr= utrs.append( [self.chrom, '.', '3utr', start_of_3utr, e, '0.000000', self.strand, -1, self.attrs] )
            if s >= start_of_3utr:
                ## Exon is fully 3utr
                xutr= utrs.append( [self.chrom, '.', '3utr', s, e, '0.000000', self.strand, -1, self.attrs] )
        self.checkCorrectCoords(utrs)
        return utrs

    def getGtf_utr_reverse(self):
        utrs= []
        if self.strand == '+':
            return utrs
        if not self.hasCDS():
            for s, e in zip(self.exonStarts, self.exonEnds):
                xutr= utrs.append( [self.chrom, '.', '3utr', s+1, e, '0.000000', self.strand, -1, self.attrs] )
            self.checkCorrectCoords(utrs)
            return utrs            
        
        end_of_3utr= self.get_stop_codon()[0][3] - 1
        start_of_5utr= self.get_start_codon()[-1][4] + 1

        for i in range(self.exonCount):
            s, e= self.exonStarts[i]+1, self.exonEnds[i]
            if e <= end_of_3utr:
                ## Exon is fully 3utr
                xutr= utrs.append( [self.chrom, '.', '3utr', s, e, '0.000000', self.strand, -1, self.attrs] )
            if s <= end_of_3utr and e >= end_of_3utr:
                ## Exon is partially a 3utr
                xutr= utrs.append( [self.chrom, '.', '3utr', s, end_of_3utr, '0.000000', self.strand, -1, self.attrs] )
            if s <= start_of_5utr and e >= start_of_5utr:
                ## Exon is partially a 5utr
                xutr= utrs.append( [self.chrom, '.', '5utr', start_of_5utr, e, '0.000000', self.strand, -1, self.attrs] )
            if s >= start_of_5utr:
                ## Exon is fully 5utr
                xutr= utrs.append( [self.chrom, '.', '5utr', s, e, '0.000000', self.strand, -1, self.attrs] )
        self.checkCorrectCoords(utrs)
        return utrs


## ---------------------------------------------------------------------
if len(sys.argv) > 1 and sys.argv[1] == '-h':
    docstring= """
Convert UCSC annotation table (genePred tables) to GTF. 

Usage:
ucscTableToGTF.py <table|stdin>

See also https://github.com/dariober/bioinformatics-cafe/tree/master/ucscTableToGTF
"""
    print docstring
    sys.exit(1)


if len(sys.argv) == 1 or sys.argv[1] == "-":
    fin= sys.stdin
else:
    fin= open(sys.argv[1])


for line in fin:
    if line.strip() == '' or line.strip().startswith('#'):
        continue
    ucsc= ucscLine(line)
    for feature in ucsc.getGtfExons():
        print '\t'.join([str(x) for x in feature])
    for feature in ucsc.getGtfIntrons():
        print '\t'.join([str(x) for x in feature])
    for feature in ucsc.get_start_codon():
        print '\t'.join([str(x) for x in feature])
    for feature in ucsc.get_stop_codon():
        print '\t'.join([str(x) for x in feature])
    for feature in ucsc.getGtfCDS():
        print '\t'.join([str(x) for x in feature])
    for feature in ucsc.getGtf_utr_forward():
        print '\t'.join([str(x) for x in feature])
    for feature in ucsc.getGtf_utr_reverse():
        print '\t'.join([str(x) for x in feature])
    for feature in ucsc.getTSS():
        print '\t'.join([str(x) for x in feature])
    for feature in ucsc.getTES():
        print '\t'.join([str(x) for x in feature])
fin.close()