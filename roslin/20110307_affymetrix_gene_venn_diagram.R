# -----------------------------------------------------------------------------
# Venn diagrams for the genes in the Affymetrix BMDM experiment
# -----------------------------------------------------------------------------

library(RODBC)
library(limma)
library(gplots)
conn<- odbcConnect(dsn= 'pgVitelleschi')

## See script 20110304_affymetrix_pig_bmdm_annotate.sql for where these data
## come from.
degenes<- sqlQuery(conn, "select * from affymetrix_bmdm_genes_de_ct", stringsAsFactors = FALSE);
degenes[1:10, ]

venndat<- degenes[, c('regulation', '0 vs 2', '0 vs 24', '0 vs 7')]
venndat[is.na(venndat)] <- FALSE
venndat[2:ncol(venndat)][venndat[, 2:ncol(venndat)] != 0] <- 1
venndat[1:10, ]

## Venn counts
vcounts_up<- vennCounts(venndat[venndat$regulation == 'upregulated', 2:ncol(venndat)])
vcounts_down<- vennCounts(venndat[venndat$regulation == 'downregulated', 2:ncol(venndat)])

names(venndat)

## Diagram for upregulated:

A = venndat[which(venndat$regulation == 'upregulated'), "0 vs 2"]
B = venndat[which(venndat$regulation == 'upregulated'), "0 vs 24"]
C = venndat[which(venndat$regulation == 'upregulated'), "0 vs 7"]

d<- list()
d$table<- table(A, B, C)
d$labels<- c("0 vs 2","0 vs 24","0 vs 7")

plot.venn.diagram(d) 

## Diagram for downregulated

A = venndat[which(venndat$regulation == 'downregulated'), "0 vs 2"]
B = venndat[which(venndat$regulation == 'downregulated'), "0 vs 24"]
C = venndat[which(venndat$regulation == 'downregulated'), "0 vs 7"]

d<- list()
d$table<- table(A, B, C)
d$labels<- c("0 vs 2","0 vs 24","0 vs 7")

plot.venn.diagram(d) 


# -----------------------------------------------------------------------------
# Script for plotting Venn diagram
# see http://tolstoy.newcastle.edu.au/R/help/03a/1115.html
# -----------------------------------------------------------------------------

venn.overlap <-
function(r, a, b, target = 0)
{
#
# calculate the overlap area for circles of radius a and b
# with centers separated by r
# target is included for the root finding code
#
        pi = acos(-1)
        if(r >= a + b) {
                return( - target)
        }
        if(r < a - b) {
                return(pi * b * b - target)
        }
        if(r < b - a) {
                return(pi * a * a - target)
        }
        s = (a + b + r)/2
        triangle.area = sqrt(s * (s - a) * (s - b) * (s - r))
        h = (2 * triangle.area)/r
        aa = 2 * atan(sqrt(((s - r) * (s - a))/(s * (s - b))))
        ab = 2 * atan(sqrt(((s - r) * (s - b))/(s * (s - a))))
        sector.area = aa * (a * a) + ab * (b * b)
        overlap = sector.area - 2 * triangle.area
        return(overlap - target)
}

plot.venn.diagram <-
function(d)
{
#
# Draw Venn diagrams with proportional overlaps
# d$table = 3 way table of overlaps
# d$labels = array of character string to use as labels
#
pi = acos(-1)
csz = 0.1
# Normalize the data
n = length(dim(d$table))
c1 = vector(length = n)
c1[1] = sum(d$table[2, , ])
c1[2] = sum(d$table[, 2, ])
c1[3] = sum(d$table[, , 2])
n1 = c1
#
c2 = matrix(nrow = n, ncol = n, 0)
c2[1, 2] = sum(d$table[2, 2, ])
c2[2, 1] = c2[1, 2]
c2[1, 3] = sum(d$table[2, , 2])
c2[3, 1] = c2[1, 3]
c2[2, 3] = sum(d$table[, 2, 2])
c2[3, 2] = c2[2, 3]
n2 = c2
#
c3 = d$table[2, 2, 2]
n3 = c3
c2 = c2/sum(c1)
c3 = c3/sum(c1)
c1 = c1/sum(c1)
n = length(c1)
# Radii are set so the area is proporitional to number of counts
pi = acos(-1)
r = sqrt(c1/pi)
r12 = uniroot(venn.overlap, interval = c(max(r[1] - r[2], r[2] - r[1],
    0) + 0.01, r[1] + r[2] - 0.01), a = r[1], b = r[2], target = c2[1, 2])$root
r13 = uniroot(venn.overlap, interval = c(max(r[1] - r[3], r[3] - r[1], 0) + 0.01,
    r[1] + r[3] - 0.01), a = r[1], b = r[3], target = c2[1, 3])$root
r23 = uniroot(venn.overlap, interval = c(max(r[2] - r[3], r[3] - r[2], 0) + 0.01,
    r[2] + r[3] - 0.01), a = r[2], b = r[3], target = c2[2, 3])$root
s = (r12 + r13 + r23)/2
x = vector()
y = vector()
x[1] = 0
y[1] = 0
x[2] = r12
y[2] = 0
angle = 2 * atan(sqrt(((s - r12) * (s - r13))/(s * (s - r13))))
x[3] = r13 * cos(angle)
y[3] = r13 * sin(angle)
xc = cos(seq(from = 0, to = 2 * pi, by = 0.01))
yc = sin(seq(from = 0, to = 2 * pi, by = 0.01))
cmx = sum(x * c1)
cmy = sum(y * c1)
x = x - cmx
y = y - cmy
rp=sqrt(x*x + y*y)
frame()
par(usr = c(-1, 1, -1, 1), pty = "s")
## box()
for(i in 1:3) {
    linecol= c('firebrick4', 'dodgerblue', 'seagreen4')
    lines(xc * r[i] + x[i], yc * r[i] + y[i], col= linecol[i], lwd= 2)
}
xl = (rp[1] + (0.7 * r[1])) * x[1]/rp[1]
yl = (rp[1] + (0.7 * r[1])) * y[1]/rp[1]
# text(xl, yl, d$labels[1])                  ## Text group name
# text(xl, yl - csz, d$table[2, 1, 1])     ## Text number 
xl = (rp[2] + (0.7 * r[2])) * x[2]/rp[2]
yl = (rp[2] + (0.7 * r[2])) * y[2]/rp[2]   ## Text number 
# text(xl, yl, d$labels[2])                  ## Text group name
# text(xl, yl - csz, d$table[1, 2, 1])
xl = (rp[3] + (0.7 * r[3])) * x[3]/rp[3]
yl = (rp[3] + (0.7 * r[3])) * y[3]/rp[3]
# text(xl, yl, d$labels[3])                  ## Text group name
# text(xl, yl - csz, d$table[1, 1, 2])
#
# text((x[1] + x[2])/2, (y[1] + y[2])/2, d$table[2, 2, 1]) ## Text number
# text((x[1] + x[3])/2, (y[1] + y[3])/2, d$table[2, 1, 2]) ## Text number
# text((x[2] + x[3])/2, (y[2] + y[3])/2, d$table[1, 2, 2]) ## Text number
#
# text(0, 0, n3)  ## Text common to all
list(r = r, x = x, y = y, dist = c(r12, r13, r23), count1 = c1, count2 =
c2, labels = d$labels)

## Group names
legend( "topright", legend= d$labels, col= linecol, bty= 'n', pch= 1,
    pt.cex= 1.5, pt.lwd= 2 )
}

#
# TRITUME
#

vennDiagram(vcounts, circle.col= c('firebrick4', 'dodgerblue', 'seagreen4'), lwd= 2)

head(degenes)

> vcounts
     0 vs 2 0 vs 24 0 vs 7 Counts
[1,]      0       0      0      0
[2,]      0       0      1    374
[3,]      0       1      0    106
[4,]      0       1      1    177
[5,]      1       0      0     27
[6,]      1       0      1     69
[7,]      1       1      0      8
[8,]      1       1      1     78
attr(,"class")
[1] "VennCounts"
:>
?venn
vcounts[-1,]


vpvenn<- venn(venndat)



