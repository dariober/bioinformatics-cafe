import re

CHROMS= ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']

wildcard_constraints:
    chrom= '|'.join([re.escape(x) for x in CHROMS])

rule all:
    input:
        'gnomad.genomes.r2.0.2.sites.hg38.simple.vcf.gz.tbi'

localrules: getLiftoverChains
rule getLiftoverChains:
    output:
        'hg19ToHg38.over.chain.gz'
    shell:
        r"""
        curl -s -O ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
        """

rule simplifyChroms:
    # From gnomAD VCF source:
    # * Keep only INFO/AF field 
    # * Rename from ENSEMBL to UCSC scheme by adding chr prefix (e.g. 22 -> chr22)
    # * Keep only biallelic sites
    output:
        vcf= temp('{chrom}.simple.vcf.gz'),
        tbi= temp('{chrom}.simple.vcf.gz.tbi'),
    params:
        cpus= 2,
        mem= 1000,
    shell:
        r"""
        bcftools annotate -O v -i 'INFO/AF > 0' -x '^INFO/AF' https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr{wildcards.chrom}.vcf.bgz \
        | awk '{{ if($1 ~ "^#") {{print $0}} else {{print "chr"$0}} }}' \
        | bcftools view -O z --min-alleles 2 --max-alleles 2 > {output.vcf}
        tabix {output.vcf}

        rm gnomad.genomes.r2.0.2.sites.chr{wildcards.chrom}.vcf.bgz.tbi
        """

rule fixVcfHeader:
    # Probably neccessary but recommended to keep downstream tools happy.
    input:
        vcf= temp('{chrom}.simple.vcf.gz'),
        tbi= temp('{chrom}.simple.vcf.gz.tbi'),
    output:
        uvcf= temp('{chrom}.fix.vcf'),
        vcf= temp('{chrom}.fix.vcf.gz'),
        tbi= temp('{chrom}.fix.vcf.gz.tbi'),
    params:
        ref= config['ref'],
        picard= config['picard'],
        cpus= 1,
        mem= 4000,
    shell:
        r"""
        java -Xmx3g -jar {params.picard} UpdateVcfSequenceDictionary I={input.vcf} O={output.uvcf} SD={params.ref}
        bgzip -c {output.uvcf} > {output.vcf}
        tabix {output.vcf}
        """

rule liftOver:
    input:
        vcf= temp('{chrom}.fix.vcf.gz'),
        tbi= temp('{chrom}.fix.vcf.gz.tbi'),
        chain= 'hg19ToHg38.over.chain.gz',
        ref= config['ref'],
    output:
        vcf= temp('{chrom}.lo.vcf.gz'),
        tbi= temp('{chrom}.lo.vcf.gz.tbi'),
        reject= temp('{chrom}.reject.vcf'),
    params:
        picard= config['picard'],
        cpus= 1,
        mem= 8000,
    shell:
        r"""
        java -Xmx6g -jar {params.picard} LiftoverVcf \
            I={input.vcf} \
            O={output.vcf} \
            CHAIN={input.chain} \
            REJECT={output.reject} \
            R={input.ref} \
            WARN_ON_MISSING_CONTIG=true
        """

rule catChroms:
    input:
        vcf= expand('{chrom}.lo.vcf.gz', chrom= CHROMS),
        tbi= expand('{chrom}.lo.vcf.gz.tbi', chrom= CHROMS),
    output:
        vcf= 'gnomad.genomes.r2.0.2.sites.hg38.simple.vcf.gz',
        tbi= 'gnomad.genomes.r2.0.2.sites.hg38.simple.vcf.gz.tbi',
    params:
        cpus= 1,
        mem= 1000,
    shell:
        r"""
        bcftools concat --allow-overlaps {input.vcf} > {output.vcf}
        tabix {output.vcf}
        """
