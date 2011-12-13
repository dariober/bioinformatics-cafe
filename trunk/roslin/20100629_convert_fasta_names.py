#! /usr/bin/python

"""

 Convert the names of the sequences in a FASTA file to be SAM compliant. I.e. not
 to contain the following characters: [ \t\n\r@=]
 These characters are replaced by '_'. A file reporting the conversion old -> new
 is generated

"""

import re

#------------------------------[ User input ]----------------------------------

fasta_filename= '/exports/work/vet_roslin_nextgen/dario/repBase/fasta/pig_repbase_20090604.fa'

output_filename= fasta_filename + '2'

dbnames_filename= fasta_filename + '.dbnames'

#------------------------------------------------------------------------------

fasta_file= open(fasta_filename)

fasta_out= open(output_filename, 'w')
dbnames_out= open(dbnames_filename, 'w')

for line in fasta_file:
    line= line.rstrip('\n')
    if line.startswith('>'):
        ori_line= line
        new_line= re.sub('[ |\t|\n|\r|@|=]', '_', line)
        db_line= [ori_line, new_line]
        fasta_out.write(new_line + '\n')
        dbnames_out.write('\t'.join(db_line) + '\n')
    else:
        fasta_out.write(line + '\n')

fasta_file.close()
fasta_out.close()
dbnames_out.close()

        


