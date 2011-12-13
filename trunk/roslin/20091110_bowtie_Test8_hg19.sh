## CAGE library Test8, filtered to have each read with no base call <20 
## aligned against human genome hg19
## Same script as 20091109_bowtie_Test8_hg19_ver2.sh but running on ebrlxserv06
## Using bowtie 0.11.3
## Allow up to 3 mismatches in the whole read (SOAP like)
## -f : Input file is fasta format
## -m 10 : Reject tags which align in more than 10 (-m 10) locations
## -k 10 : For each tag, report up to 10 locations (i.e all of those allowed by -m)

/nfs_netapp/dberaldi/bowtie/bowtie-0.11.3/bowtie h_sapiens_37_asm -q /nfs_netapp/dberaldi/Documents/fastq/s_6_test8_qs20.fastq -m 10 -k 10 -v 3 /nfs_netapp/dberaldi/Documents/Tritume/bwt_Test8_hg19_v3.map