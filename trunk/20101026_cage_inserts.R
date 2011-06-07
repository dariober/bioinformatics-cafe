no_mism<- c(12, 6, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0)
position<- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27)

barplot(no_mism, names= position, 
    xlab= 'Position on CAGE tag (5\' -> 3\')',
    ylab= 'No mismatches', main= 'Where mismatches occur, (n= 29 mismatches)')
savePlot('M:/Documents/LabBook/LabBook_Figures/20101027_mismatch_pos_cage.emf', 'emf')

bases<- read.table('D:/Tritume/fasta_to_tabular.txt', header=T, sep='\t')
bases[1:10,]

base_by_pos<- aggregate(bases$base, by= list(bases$pos, bases$base), length)
names(base_by_pos)<- c('pos', 'base', 'count')

basew<- reshape(data.frame(base_by_pos), v.names="count", idvar="pos",
                timevar="base", direction="wide")
rownames(basew)<- basew[,1]
basew<- basew[, -1]

barplot(t(basew), xlab= 'Position on CAGE tag (5\' -> 3\')', 
    ylab= 'Base count', legend.text = c('A', 'C', 'G', 'T'), 
    args.legend= list(x= 34.5, bty= 'n'), 
    col= c('green', 'steelblue1', 'grey29', 'firebrick3'))
