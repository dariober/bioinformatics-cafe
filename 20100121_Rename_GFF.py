# -----------------[Users input]-----------------------------------------------

## Map current chromosomes names (left) to new ones (right):

map_chr_names= {
'MT' :'MT dna:chromosome chromosome:Sscrofa9:MT:1:16613:1', 
'18' :'18 dna:chromosome chromosome:Sscrofa9:18:1:54314914:1', 
'12' :'12 dna:chromosome chromosome:Sscrofa9:12:1:57436344:1', 
'17' :'17 dna:chromosome chromosome:Sscrofa9:17:1:64400339:1', 
'10' :'10 dna:chromosome chromosome:Sscrofa9:10:1:66741929:1', 
'16' :'16 dna:chromosome chromosome:Sscrofa9:16:1:77440658:1', 
'11' :'11 dna:chromosome chromosome:Sscrofa9:11:1:79819395:1', 
'5' :'5 dna:chromosome chromosome:Sscrofa9:5:1:100521970:1', 
'8' :'8 dna:chromosome chromosome:Sscrofa9:8:1:119990671:1', 
'6' :'6 dna:chromosome chromosome:Sscrofa9:6:1:123310171:1', 
'3' :'3 dna:chromosome chromosome:Sscrofa9:3:1:123604780:1', 
'X' :'X dna:chromosome chromosome:Sscrofa9:X:1:125876292:1', 
'9' :'9 dna:chromosome chromosome:Sscrofa9:9:1:132473591:1', 
'15' :'15 dna:chromosome chromosome:Sscrofa9:15:1:134546103:1', 
'4' :'4 dna:chromosome chromosome:Sscrofa9:4:1:136259946:1', 
'7' :'7 dna:chromosome chromosome:Sscrofa9:7:1:136414062:1', 
'2' :'2 dna:chromosome chromosome:Sscrofa9:2:1:140138492:1', 
'13' :'13 dna:chromosome chromosome:Sscrofa9:13:1:145240301:1', 
'14' :'14 dna:chromosome chromosome:Sscrofa9:14:1:148515138:1', 
'1' :'1 dna:chromosome chromosome:Sscrofa9:1:1:295534705:1', 
'human_rDNA' :'gi|555853|gb|U13369.1|HSU13369 Human ribosomal DNA complete repeating unit'}

## Input gff3 file
infile= 'C:/Tritume/Sus_scrofa.Sscrofa9.56.gff3'

## Output is the same as input with 'ensnames appended

#------------------------------------------------------------------------------

gff_ss_handle= open(infile)
gff_ss_out= open(infile + '.ensnames', 'w')

for line in gff_ss_handle:
    if line.startswith('##'):
        gff_ss_out.write(header_line.rstrip('\n') + '; Dario: Chromosome names changed to ensembl long format\n')
        continue
    xlist= line.split('\t')
    chr_name= xlist[0]
    chr_name_out= map_chr_names[chr_name]
    xlist[0] = chr_name_out
    gff_ss_out.write('\t'.join(xlist))
    
gff_ss_handle.close()
gff_ss_out.close()