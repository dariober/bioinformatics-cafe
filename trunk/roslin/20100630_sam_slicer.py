import time

filename_sam= 'D:/Tritume/accepted_hits.sam'
filename_slice= filename_sam + '.slice'

contig= '4'

sam= open(filename_sam)
slice= open(filename_slice, 'w')

t0=time.time()
counter= 0
print('\nReading SAM file...\n')
for line in sam:
    line= line.split('\t')
    if line[2] == contig:
        line= '\t'.join(line)
        slice.write(line)
    counter += 1
    if counter % 1000000 == 0:
        print('Read SAM line number: ' + str(counter))
t1= time.time()
print('\nFinished reading SAM file in '+ str(round(t1-t0, 2)) + ' sec.\n')
print(str(counter) + ' lines read.\n')

sam.close()
slice.close()
