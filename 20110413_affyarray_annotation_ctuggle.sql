select read_table($$ file: 'D:/Tritume/toChris_Annotation_UniqID.txt', table: 'affymetrix.annotation_ctuggle', overwrite:True,
apply: """
def apply(line):
    line= line.split('\t')
    gene= line[1]
    gene= gene.split('.')
    gene_name= '.'.join(gene[:-1])
    if gene_name == 'UNANN-':
        gene_name= ''
    gene_seq= gene[-1]
    line= line + [gene_name] + [gene_seq]
    line= '\t'.join(line)
    return(line)
""", skip: 1, header: ['probe_set_id', 'data_id', 'gene_name', 'gene_seq']  $$)

comment on table affymetrix.annotation_ctuggle is 'Annotation for the porcine 24K array from Chris Tuggle. see 20110413_affyarray_annotation_ctuggle.sql. Original data is toChris_Annotation_UniqID.xlsx'

