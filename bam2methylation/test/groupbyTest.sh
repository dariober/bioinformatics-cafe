cd ~/svn_git/bioinformatics-cafe/trunk/bam2methylation

cat test_data/groupby.bed | bedtools groupby -g 1,2,3 -c 5,6,7,8 -o sum,sum,distinct,sum > exp.tmp
cat test_data/groupby.bed | python groupby.py > obs.tmp
dx=`diff exp.tmp obs.tmp`

if [[ "$dx" != "" ]]
then
    echo "** FAILED **"
    diff exp.tmp obs.tmp
fi
rm exp.tmp obs.tmp

echo "CAN GROUP ONE LINE ONLY"
head -n 1 test_data/groupby.bed | bedtools groupby -g 1,2,3 -c 5,6,7,8 -o sum,sum,distinct,sum > exp.tmp
head -n 1 test_data/groupby.bed | python groupby.py > obs.tmp
dx=`diff exp.tmp obs.tmp`

if [[ "$dx" != "" ]]
then
    echo "** FAILED **"
    diff exp.tmp obs.tmp
fi
rm exp.tmp obs.tmp

echo "CAN GROUP LAST LINE SINGLE"
grep 'chr1' test_data/groupby.bed | bedtools groupby -g 1,2,3 -c 5,6,7,8 -o sum,sum,distinct,sum > exp.tmp
grep 'chr1' test_data/groupby.bed | python groupby.py > obs.tmp
dx=`diff exp.tmp obs.tmp`

if [[ "$dx" != "" ]]
then
    echo "** FAILED **"
    diff exp.tmp obs.tmp
fi
rm exp.tmp obs.tmp

echo "CAN GROUP ONLY SINGLE LINES"
grep '1st' test_data/groupby.bed | bedtools groupby -g 1,2,3 -c 5,6,7,8 -o sum,sum,distinct,sum > exp.tmp
grep '1st' test_data/groupby.bed | python groupby.py > obs.tmp
dx=`diff exp.tmp obs.tmp`

if [[ "$dx" != "" ]]
then
    echo "** FAILED **"
    diff exp.tmp obs.tmp
fi
rm exp.tmp obs.tmp
