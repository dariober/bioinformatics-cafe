#!/usr/bin/python

"""
Read a pileup file of aligned CAGE tags to produce tss clusters

Should look like this:

rname   pos     ref_base    read_count      strand
MT      28      G           8               +
MT      29      C           8               +
MT      30      A           8               +
MT      31      A           8               +
MT      32      A           8               +
MT      33      C           8               +


"""
import sys
import os

# -------------------------------[ Input/Output ]------------------------------

os.chdir('/exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment')

## Input/Output files:
pile_in= 'cage_050810_bwa_20101218.clean.single.pileup'
pile_out= 'cage_050810_bwa_20101218.clean.single.tss.pileup'

## Edit lines about header according to whether input has header or not

# -----------------------------------------------------------------------------

def ending():
    """ Do these things once the last line is found  """
    print(str(n) + ' positions in pileup file.')
##    print(str(n_tc) + ' TSS clusters found.')
    pileup_in.close()
    pileup_out.close()
    sys.exit()

def new_tag():
    """ Generate a new TSS tag like this:
    SSTSC_1_F_000000123
    """
    ## Convert +/- to F/R
    if line[4] == '+':
        strand= 'F'
    elif line[4] == '-':
        strand= 'R'
    else:
        sys.exit('Unrecognized character for strand: ' + strand)
    tss_tag= 'SSTSC_' + line[0] + '_' + strand + '_' + format(int(line[1]), '09d')
    return(tss_tag)

# -----------------------------------------------------------------------------

pileup_in= open(pile_in)
pileup_out= open(pile_out, 'w')

## Write out header:
h= pileup_in.readline().rstrip('\n')
pileup_out.write(h + '\ttss_cluster_id\n')

n= 0     ## counter of positions in pileup
n_tc= 0  ## counter of TSS cluster produced

## First TSS cluster ID: In the form like: 'SSTSC'_1_F_00000123
line= pileup_in.readline().rstrip('\n')
line= line.split('\t')
tss_tag= new_tag()
n_tc += 1
while True:
    if line == '':
        ending()
    p= (line[0], int(line[1]), line[4])
    line.append(tss_tag)
    pileup_out.write('\t'.join(line) + '\n')
    n += 1
    line= pileup_in.readline().rstrip('\n')
    line= line.split('\t')
##    cur_pos= (line[0], int(line[1])-1, line[4])
    cur_pos= (line[0], int(line[1]), line[4])
    
    ## Condition to keep adjacent positions in pileup under the same TSS cluster
    ##    |----Same rname----|         |-- Same position or immediately the next --|          |----Same strand----|
    while  cur_pos[0] == p[0]    and   (cur_pos[1] == p[1] or cur_pos[1] == p[1]+1)     and    cur_pos[2] == p[2]:
        p= (line[0], int(line[1]), line[4])
        line.append(tss_tag)
        pileup_out.write('\t'.join(line) + '\n')
        line= pileup_in.readline().rstrip('\n')
        n += 1
        if line == '':
            ending()
        line= line.split('\t')
        cur_pos= (line[0], int(line[1]), line[4])
    tss_tag= new_tag()
    n_tc += 1
