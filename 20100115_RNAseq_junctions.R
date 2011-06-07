library(reshape)

lps_junx<- read.table("C:/Tritume/LPS_junctions.bed", skip= 1, sep= "\t", 
  header= F)

bed_names<- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
  "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")

names(lps_junx)<- bed_names

junx_chr<- aggregate(lps_junx$name, list(lps_junx$chrom), length)
names(junx_chr)<- c("chrom", "no_junx")
junx_chr$chrom<- as.character(junx_chr$chrom)
junx_chr$chrom[junx_chr$chrom == "X"]<- "99"
junx_chr$chrom<- as.numeric(as.character(junx_chr$chrom))

junx_chr<- junx_chr[order(junx_chr$chrom), ]


  
junx_size_genome<- lps_junx$chromEnd - lps_junx$chromStart  

junx_size_blocks<- colsplit(lps_junx$blockSizes, split= ",",
  names= c("blockSize1", "blockSize2"))

intron_size<- (lps_junx$chromEnd - lps_junx$chromStart) - 
  rowSums(junx_size_blocks)

par(mfrow= c(3,2), cex= 0.75)
  hist(lps_junx$score, main= "No. of supporting alignments/junction",
    xlab= "No. alignments\n(all)")
  hist(lps_junx$score[lps_junx$score <= 10], breaks= seq(0, 10), 
    main= "", xlab= "No. alignments\n(<= 10)")
  hist(intron_size, main= "Intron size", xlab= "Size (bp)\n(all)")    
  hist(intron_size[intron_size<= 1000], main= "", 
    xlab= "Size (bp)\n(<= 1000 bp)")    
  hist(rowSums(junx_size_blocks), 
    main= "Spliced junction span\n(length of exon_end + exon_start)",
    xlab= "Size (bp)")
  barplot(junx_chr$no_junx, names= junx_chr$chrom, xlab= "Chromosome", 
    ylab= "No. junctions", main= "RNAseq_LPS\nJunctions per chromosome",
    sub= "TopHat 14/01/2010")
  mtext(text= paste("Total:", nrow(lps_junx), "junctions"), side= 3, line= -1.5,
    adj= 1, cex= 0.75, font= 2)

 

summary(intron_size)
# ------------------------------------[ Tritume ]------------------------------
max(lps_junx$blockCount)  
aggregate(lps_junx$name, list(lps_junx$strand), length)
lps_junx$score[1:10]
lps_junx[1:10,]