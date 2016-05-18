## expected. Get genePredToGtf from ucsc utilities
curl -o - -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vegaGene.txt.gz | gunzip -c \
| tail -n+2 | cut -f 2- \
| genePredToGtf -utr file stdin stdout \
| sort -k1,1 -k4,4n -k5,5n > genePredToGtf.tmp.gtf

## Observed
curl -o - -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vegaGene.txt.gz | gunzip -c \
| ./ucscTableToGTF.py \
| sort -k1,1 -k4,4n -k5,5n > ucscTableToGTF.tmp.gtf

grep 'start_codon' genePredToGtf.tmp.gtf | cut -f 1,3,4,5,7,8 > exp.tmp.txt
grep 'start_codon' ucscTableToGTF.tmp.gtf | cut -f 1,3,4,5,7,8 > obs.tmp.txt
diff exp.tmp.txt obs.tmp.txt


grep 'CDS' genePredToGtf.tmp.gtf | cut -f 1,3,4,5,7,8 > exp.tmp.txt
grep 'CDS' ucscTableToGTF.tmp.gtf | cut -f 1,3,4,5,7,8 > obs.tmp.txt
diff exp.tmp.txt obs.tmp.txt

## NOT USED ANY MORE
## =================

## Collect all genes so we test them one by one
transcripts_pos=`grep -w '+' test_data/hg19.chr19.vegaGene.gtf | sed -e 's/.*transcript_id "//' -e 's/".*//' | sort | uniq`
transcripts_neg=`grep -w '-' test_data/hg19.chr19.vegaGene.gtf | sed -e 's/.*transcript_id "//' -e 's/".*//' | sort | uniq`

## echo $transcripts | wc                                ## 1532
## grep -v '#' test_data/hg19.chr19.vegaGene.txt | wc -l ## 1532


echo "EXONS +"

for trx in $transcripts_pos
do
    grep -w ${trx} test_data/hg19.chr19.vegaGene.gtf | grep -w 'exon' | cut -f 1,3-8 > ${trx}.exp.tmp.gtf
    grep -w ${trx} test_data/hg19.chr19.vegaGene.txt | python ucscTableToGTF.py - | grep -w 'exon' | cut -f 1,3-8 > ${trx}.obs.tmp.gtf
    xdiff=`diff ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf`
    if [ "$xdiff" != "" ]
        then
        echo "ERROR in $trx"
        cat ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf | column -t
        diff ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf
        break
    fi
    rm ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf
done

echo "EXONS -"

for trx in $transcripts_neg
do
    grep -w ${trx} test_data/hg19.chr19.vegaGene.gtf | grep -w 'exon' | cut -f 1,3-8 > ${trx}.exp.tmp.gtf
    grep -w ${trx} test_data/hg19.chr19.vegaGene.txt | python ucscTableToGTF.py - | grep -w 'exon' | cut -f 1,3-8 > ${trx}.obs.tmp.gtf
    xdiff=`diff ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf`
    if [ "$xdiff" != "" ]
        then
        echo "ERROR in $trx"
        cat ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf | column -t
        diff ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf
        break
    fi
    rm ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf
done

echo "STOP_CODON +"

for trx in $transcripts_pos
do
    echo $trx
    grep -w ${trx} test_data/hg19.chr19.vegaGene.gtf | grep -w 'stop_codon' | cut -f 1,3-8 > ${trx}.exp.tmp.gtf
    grep -w ${trx} test_data/hg19.chr19.vegaGene.txt | python ucscTableToGTF.py - | grep -w 'stop_codon' | cut -f 1,3-8 > ${trx}.obs.tmp.gtf
    xdiff=`diff ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf`
    if [ "$xdiff" != "" ]
        then
        echo $trx
        echo "Expected"
        cat ${trx}.exp.tmp.gtf
        echo "Observed"
        cat ${trx}.obs.tmp.gtf
        diff ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf
        break
    fi
    rm ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf
done

echo "STOP_CODON -"

for trx in $transcripts_neg
do
    echo $trx
    grep -w ${trx} test_data/hg19.chr19.vegaGene.gtf | grep -w 'stop_codon' | cut -f 1,3-8 > ${trx}.exp.tmp.gtf
    grep -w ${trx} test_data/hg19.chr19.vegaGene.txt | python ucscTableToGTF.py - | grep -w 'stop_codon' | cut -f 1,3-8 > ${trx}.obs.tmp.gtf
    xdiff=`diff ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf`
    if [ "$xdiff" != "" ]
        then
        echo $trx
        echo "Expected"
        cat ${trx}.exp.tmp.gtf
        echo "Observed"
        cat ${trx}.obs.tmp.gtf
        diff ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf
        break
    fi
    rm ${trx}.exp.tmp.gtf ${trx}.obs.tmp.gtf
done

echo "START_CODON +"

for trx in $transcripts_pos
do
    echo $trx
    exp=`grep -w ${trx} test_data/hg19.chr19.vegaGene.gtf | grep -w 'start_codon' | cut -f 1,3-8`
    obs=`grep -w ${trx} test_data/hg19.chr19.vegaGene.txt | python ucscTableToGTF.py - | grep -w 'start_codon' | cut -f 1,3-8`
    if [ "$exp" != "$obs" ]
        then
        echo $trx
        echo "Expected"
        echo $exp
        echo "Observed"
        echo $obs
        break
    fi
done


echo "START_CODON -"

for trx in $transcripts_neg
do
    echo $trx
    exp=`grep -w ${trx} test_data/hg19.chr19.vegaGene.gtf | grep -w 'start_codon' | cut -f 1,3-8`
    obs=`grep -w ${trx} test_data/hg19.chr19.vegaGene.txt | python ucscTableToGTF.py - | grep -w 'start_codon' | cut -f 1,3-8`
    if [ "$exp" != "$obs" ]
        then
        echo $trx
        echo "Expected"
        echo $exp
        echo "Observed"
        echo $obs
        break
    fi
done

echo "CDS + "

for trx in $transcripts_pos
do
    echo $trx
    exp=`grep -w ${trx} test_data/hg19.chr19.vegaGene.gtf | grep -w 'CDS' | cut -f 1,3-8`
    obs=`grep -w ${trx} test_data/hg19.chr19.vegaGene.txt | python ucscTableToGTF.py - | grep -w 'CDS' | cut -f 1,3-8`
    if [ "$exp" != "$obs" ]
        then
        echo $trx
        echo "Expected"
        echo $exp
        echo "Observed"
        echo $obs
        break
    fi
done


echo "CDS -"

for trx in $transcripts_neg
do
    echo $trx
    exp=`grep -w ${trx} test_data/hg19.chr19.vegaGene.gtf | grep -w 'CDS' | cut -f 1,3-8`
    obs=`grep -w ${trx} test_data/hg19.chr19.vegaGene.txt | python ucscTableToGTF.py - | grep -w 'CDS' | cut -f 1,3-8`

    if [ "$exp" != "$obs" ]
        then
        echo $trx
        echo "Expected"
        echo $exp
        echo "Observed"
        echo $obs
        break
    fi
done
