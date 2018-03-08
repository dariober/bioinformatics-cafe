Effect of GATK base quality recalibration
=========================================

We want to see the effect of BQSR on a hypothetical hypermutated genome

* Copy reference sequence chr18 and add a mismatch every 1000 bp

* With a tiling approach extract reads to make 80x coverage of the entire chrom

* Extract reads from the original reference to make 20x coverage

* Combine the two sets of reads. Now you have a library sequenced at 100x depth 
  from a hypermutated sample

* Map to reference chr18

* Apply BQSR

* See the effect of BQSR

```
python -c "
import pyfaidx
import random
import textwrap

mut_dict= {'A': 'CTG', 'C': 'ATG', 'G': 'ACT', 'T': 'ACG'}

def random_baseq(n):
    q= ''.join([random.choice(["'", '0', '7', '<', 'B', 'F', 'I']) for _ in range(0, n)])
    return q

ref= pyfaidx.Fasta('/scratch/dberaldi/ref/GRCh38/GRCh38.fa')
chr18= ref['chr18'][0:].seq.strip('N').upper()
mut= []
for i in range(0, len(chr18)):
    b= chr18[i]
    if b != 'N' and i % 1000 == 0:
       m= random.randint(0, 2)
       mut.append(mut_dict[b][m])
       print(i)
    else:
        mut.append(b)

mut18= ''.join(mut)

fout= open('mut.fq', 'w')

rstart= 0
rlen= 101
rend= rlen
rnum= 1
while rend < len(mut18):
    rseq= mut18[rstart:rend]
    rstart= rend
    rend= rstart + rlen
    for i in range(0, 24):
        fq= ['@r' + str(rnum), rseq, '+', 'I' * rlen]
        rnum += 1
        fout.write('\n'.join(fq) + '\n')
    if rnum % 10000 == 0:
        print(rnum)
fout.close()
"
```

Now do the mapping:

```
samtools faidx /scratch/dberaldi/ref/GRCh38/GRCh38.fa chr18 > chr18.fa
java -jar ~/applications/picard/picard.jar CreateSequenceDictionary R=chr18.fa

bwa index chr18.fa
bwa mem -R r'@RG\tID:mut\tLB:mut\tSM:mut\tPL:ILLUMINA\tPU:NA' -t 12 chr18.fa mut.fq | samtools sort -@ 8 > mut.bam

#java -jar ~/applications/picard/picard.jar AddOrReplaceReadGroups I=mut.bam O=mut2.bam ID=mut LB=mut PL=ILLUMINA SM=mut PU=NA
#mv mut2.bam mut.bam
samtools index mut.bam
```

And recalibrate

```
java -jar ~/applications/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
           -T BaseRecalibrator \
           -R chr18.fa \
           --knownSites /scratch/dberaldi/ref/GRCh38/CosmicMuts.gatk.vcf.gz \
           -I mut.bam \
           -o mut.bqsr

java -jar ~/applications/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
     -T PrintReads \
     -R chr18.fa \
     -I mut.bam \
     -BQSR mut.bqsr \
     -o mut.bqsr.bam
```
