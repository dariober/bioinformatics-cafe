fastaf= open('M:/Documents/LabBook/ABI_sequences/20101020/cage_inserts.fa', 'r')
fasta= fastaf.readlines()
outf= open('D:/Tritume/fasta_to_tabular.txt', 'w')

outf.write('base\t' + 'tag_name\t' + 'pos\n')

for i in fasta:
    if i.startswith('>'):
        name= i.lstrip('>')
        name= name.rstrip('\n')
        continue
    n= 1
    seq= i.rstrip('\n')
    for b in seq:
        outf.write(b + '\t' + name + '\t' + str(n) + '\n')
        n += 1

outf.close()
fastaf.close()