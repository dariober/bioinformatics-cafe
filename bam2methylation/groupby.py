#Python equivalent of bedtools command:
#bedtools groupby -g 1,2,3 -c 5,6,7,8 -o sum,sum,distinct,sum
#USAGE: No argument taken. Read from stdin, write to stdout

import sys

def groupToString(group):
    return '\t'.join([group[0], group[1], group[2], str(group[4]), str(group[5]), '\t'.join(list(group[6])), str(group[7])])

group= None
prev= None
for curr in sys.stdin:
    curr= curr.strip().split('\t')
    if len(curr) < 8:
        sys.stderr.write("Number of fields < 8 in line " + curr)
        print ''
        sys.exit(1)
    curr[4]= int(curr[4])
    curr[5]= int(curr[5])
    curr[7]= int(curr[7])
    if not group:
        group= [x for x in curr[0:8]]
        group[6]= set(curr[6])
    if not prev:
        pass
    elif curr[0:3] == prev[0:3]:
        group[4] += curr[4]
        group[5] += curr[5]
        group[6].add(curr[6])
        group[7] += curr[7]
    elif curr[0:3] != prev[0:3]:
        print groupToString(group)
        group= [x for x in curr[0:8]]
        group[6]= set(curr[6])
    else:
        sys.stderr.write('Unexpected state'); print ''
        sys.exit(1)
    prev= curr

print groupToString(group)
    
sys.exit()


