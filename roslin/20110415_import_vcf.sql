/*
    Import vcf file produced from samtools mpileup | bcftools from RNAseq am and bmdm.
    See Eddie /exports/work/vet_roslin_nextgen/dario/samtools/output/20110414_rnaseq_mpileup for pileups and filetring.
    See /exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_xxx for alignments.

    This python script was used to prepare the SQL statement (manually edited to have DROP TABLE and append:False)

fout= open('D:/Tritume/import_vcf.sql', 'w')
for file_name in ('var.flt.am_cat.vcf',  'var.flt.am_ctrl.vcf', 'var.flt.am_lps.vcf',  'var.flt.bmdm_cat.vcf',  'var.flt.bmdm_ctrl.vcf',  'var.flt.bmdm_lps.vcf'):
    read_table= """
select read_table($$ file: 'D:/Tritume/%s', table: 'vcf_variants.vcf', comment_char: '#', append:True,
    apply:\"\"\"def apply(line):
                        if line.startswith('#'):
                            return(line)
                        line= '%s' + '\t' + line
                        return(line)\"\"\",
    header: ["file", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "genotype"] $$);\n""" %(file_name, file_name)
    fout.write(read_table)
fout.close()
   
*/

-- drop table vcf_variants.vcf;
select read_table($$ file: 'D:/Tritume/var.flt.am_cat.vcf', table: 'vcf_variants.vcf', comment_char: '#', append:False,
    apply:"""def apply(line):
                        if line.startswith('#'):
                            return(line)
                        line= 'var.flt.am_cat.vcf' + '	' + line
                        return(line)""",
    header: ["file", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "am", "bmdm"] $$);

select read_table($$ file: 'D:/Tritume/var.flt.am_ctrl.vcf', table: 'vcf_variants.vcf', comment_char: '#', append:True,
    apply:"""def apply(line):
                        if line.startswith('#'):
                            return(line)
                        line= 'var.flt.am_ctrl.vcf' + '	' + line
                        return(line)""",
    header: ["file", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "am", "bmdm"] $$);

select read_table($$ file: 'D:/Tritume/var.flt.am_lps.vcf', table: 'vcf_variants.vcf', comment_char: '#', append:True,
    apply:"""def apply(line):
                        if line.startswith('#'):
                            return(line)
                        line= 'var.flt.am_lps.vcf' + '	' + line
                        return(line)""",
    header: ["file", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "am", "bmdm"] $$);

select read_table($$ file: 'D:/Tritume/var.flt.bmdm_cat.vcf', table: 'vcf_variants.vcf', comment_char: '#', append:True,
    apply:"""def apply(line):
                        if line.startswith('#'):
                            return(line)
                        line= 'var.flt.bmdm_cat.vcf' + '	' + line
                        return(line)""",
    header: ["file", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "am", "bmdm"] $$);

select read_table($$ file: 'D:/Tritume/var.flt.bmdm_ctrl.vcf', table: 'vcf_variants.vcf', comment_char: '#', append:True,
    apply:"""def apply(line):
                        if line.startswith('#'):
                            return(line)
                        line= 'var.flt.bmdm_ctrl.vcf' + '	' + line
                        return(line)""",
    header: ["file", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "am", "bmdm"] $$);

select read_table($$ file: 'D:/Tritume/var.flt.bmdm_lps.vcf', table: 'vcf_variants.vcf', comment_char: '#', append:True,
    apply:"""def apply(line):
                        if line.startswith('#'):
                            return(line)
                        line= 'var.flt.bmdm_lps.vcf' + '	' + line
                        return(line)""",
    header: ["file", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "am", "bmdm"] $$);

comment on table vcf is 'vcf files produced by samtool mpileup | bcftools see /exports/work/vet_roslin_nextgen/dario/samtools/output/20110414_rnaseq_mpileup.
  var.flt.bmdm_cat.vcf   -> both bmdm files concatenated and send to mpileup
  var.flt.am_cat.vcf     -> both am files concatenated and send to mpileup
  var.flt.bmdm_lps.vcf   -> Individual alignments
  var.flt.am_ctrl.vcf 
  var.flt.am_lps.vcf
  var.flt.bmdm_ctrl.vcf 
Alignments comes from 20110202_rnaseq_am_ctrl/accepted_hits.bam & 20110202_rnaseq_am_ctrl/accepted_hits.bam; 20110202_rnaseq_bmdm_ctrl/accepted_hits.bam & 20110202_rnaseq_bmdm_ctrl/accepted_hits.bam';


