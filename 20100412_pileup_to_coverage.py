#!/usr/bin/python

"""
Returns the number of positions (count_of_positions) covered by coverage_depth reads.

Input file is in pileup format where the 8th field (tab separated) is the number of
reads covering that position

-----------------------------[ User's input ]--------------------------------- """

pileup_infile= '/exports/work/vet_roslin_nextgen/dario/samtools/output/20100409_RNAseq_CTRL_tophat/20100409_RNAseq_CTRL_sscrofa9.56.pileup'

cov_outfile= pileup_infile + '.cov'

# -----------------------------------------------------------------------------

pileup= open(pileup_infile, 'r')
cov_file= open(cov_outfile, 'w')

cov_file.write('coverage_depth\tcount_of_positions\n')

cov_dict={} 
counter= 0
for line in pileup:
    line= line.split('\t')
    cov= int(line[7]) ## len(format_allele(line))
    cov_dict[cov]= cov_dict.get(cov, 0) + 1
    counter += 1
    if counter % 5000000 == 0:
        print(str(counter) + ' rows processed')  
##        break

print(str(counter) + ' positions found in ' + pileup_infile)
print('Writing out data...')

cov_list= list(cov_dict.items())
cov_list= [list(x) for x in cov_list]
cov_list.sort()
for count in cov_list:
    out_line= '\t'.join([str(x) for x in count]) + '\n'
    cov_file.write(out_line)

pileup.close()
cov_file.close()

    

    
