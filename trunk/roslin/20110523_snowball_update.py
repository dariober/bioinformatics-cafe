#
# Prepare SQL statement to update blastn_nonpsr_refseqrna_species
# This script linked to 20110504_blast_snowball-2.0.sql
#

spec_order= ['Homo sapiens', 'Sus scrofa', 'Bos taurus', 'Pan troglodytes',
             'Mus musculus', 'Canis lupus familiaris', 'Pongo abelii',
             'Equus caballus', 'Rattus norvegicus', 'Macaca mulatta']

print("-- Not LOC annotation")
for s in spec_order:
    template= """update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = '%s' and symbol not like 'LOC%%' and evalue < 1e-9 and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);\n""" %(s)
    print(template)

print("""-- All remaining not LOC genes:
update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name like '%' and symbol not like 'LOC%' and bitscore_rank = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);\n""")

print("\n-- Use LOC annotation, whatever the evalue:")
for s in spec_order:
    template= """update blastn_nonpsr_refseqrna_species set approved_hit= True
where species_name = '%s' and bitscore_rank_species = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);\n""" %(s)
    print(template)

print("""-- Use everything else:
update blastn_nonpsr_refseqrna_species set approved_hit= True
where bitscore_rank = 1 and
qseqid not in (select distinct qseqid from blastn_nonpsr_refseqrna_species where approved_hit is true);\n""")

