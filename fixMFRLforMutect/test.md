
Fails with gatk 4.0.1.0:

```
gatk FilterMutectCalls --variant test_data/in.vcf --output out.vcf.gz
rm out.vcf.gz out.vcf.gz.tbi
```

Fix it

```
./fixMFRLforMutect.py test_data/in.vcf > fix.vcf

grep '123400000,471' fix.vcf
```

Success

```
gatk FilterMutectCalls --variant fix.vcf --output out.vcf.gz
rm out.vcf.gz out.vcf.gz.tbi
rm fix.vcf
```

