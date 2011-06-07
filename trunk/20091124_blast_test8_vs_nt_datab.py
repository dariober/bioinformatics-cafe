inblast= open('M:/Tritume/Blast_Test8_qs20_not_aligned_against_nt_db.txt', 'r')
best_hits= open('M:/Tritume/best_hits.txt', 'w')

while 1<2:
    line= inblast.readline()
    if 'Sequences producing' in line:
        inblast.readline()
        best_line= inblast.readline()
 ##       best_line= best_line[(best_line.find('|')+1) : best_line.rfind('|')] + '\n'
        best_hits.write(best_line)
    if len(line) == 0:
        break
    
inblast.close()
best_hits.close()