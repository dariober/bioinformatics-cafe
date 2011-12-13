#!/usr/bin/python

sam_all= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101222_blackface_dge/LN_all.single.sam')
tag_file= open('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101222_blackface_dge/blackface_tag_count.csv', 'w')

## LN114   HWI-EAS222_30JTBAAXX:7:1:100:1771       16      12      15227703        23

group_dict= {}
n= 0
for line in sam_all:
    line= line.split('\t')
    if int(line[5]) <= 15:
        continue
    tag= '\t'.join([line[0], line[3], line[4], line[2]])
    group_dict[tag]= group_dict.get(tag, 0)+1
    n += 1
    if n % 1000000 == 0:
        print(str(n) + ' lines to dictionary')
print(str(n) + ' lines to dictionary')
print('Writing out dictionary...')
sam_all.close()
tag_file.write('library_id\trname\tpos\tstrand\ttag_count\n')
for k in group_dict:
    tag_file.write(k + '\t' + str(group_dict[k]) + '\n')
tag_file.close()