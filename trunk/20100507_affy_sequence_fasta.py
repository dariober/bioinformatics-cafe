"""
Convert the annotation file for the porcine chip na30 (affymetrix)
to FASTA formt

Input file looks like this:

Probe Set ID	probe x	probe y	probe interrogation position	probe sequence	target strandedness		
RPTR-Ssc-J01636-4_at	666	691	118	TTTAATCACTCGCATCCATCAGAAG	Antisense		
RPTR-Ssc-J01636-4_at	511	587	213	TGTCTATTTCTCTTACGGTTCCAAC	Antisense		
RPTR-Ssc-J01636-4_at	622	447	230	GTTCCAACATCCATATAGGCCGCAA	Antisense		
RPTR-Ssc-J01636-4_at	180	587	309	TGATAACGTACTGATTGCACCCAAC	Antisense		
RPTR-Ssc-J01636-4_at	33	547	352	GGACACCCTGTACACCATGAATTGA	Antisense		

"""

file_in= 'U:/Documents/affimetrix/Porcine.probe_tab'

fasta_in= open(file_in, 'r')
file_out= 'U:/Documents/affimetrix/Porcine_na30.fa'
fasta_out= open(file_out, 'w')

fasta_in.readline()

for line in fasta_in:
    line= line.split('\t')
    fa_name= ['>', 'probe_set_id:', line[0], '|', 'probe_id:', line[0], '_', line[1], '_', line[2], '\n']
    fa_name= ''.join(fa_name)
    fasta_out.write(fa_name)
    fa_seq= line[4] + '\n'
    fasta_out.write(fa_seq)

fasta_in.close()
fasta_out.close()




