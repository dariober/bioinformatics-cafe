bedDist<- function(start1, end1, start2, end2, strand= NA){
    "
    DESCRIPTION
    Compute distance between bed features *1 and *2.
    - If strand of feature 1 is + and feature 2 has lower coords, the distance is -ve
    - If strand of feature 1 is - and feature 2 has larger coords, the distance is -ve
    - If 1 and 2 overlap distance is 0
    
    ARGS
    start1, end1, start2, end2: Vectors of start and end positions
    strand: Vector of strand of feature 1 (default +). Needs to be + or -
    
    EXAMPLES:
    bedDist(start1= 1, end1= 10, start2= 15, end2= 30) ==> 5
    bedDist(start1= 15, end1= 30, start2= 1, end2= 10) ==> -5
    bedDist(start1= 1, end1= 10, start2= 15, end2= 30, strand= '-') ==> -5
    bedDist(start1= 1, start2= 5, end1= 10, end2= 30) ==> 0
    bedDist(start2= 1, start1= 5, end2= 10, end1= 30) ==> 0"
    
    if (is.na(strand)[1]){
        strand<- rep('+', length(start1))
    }
    arglen<- c(length(start1), length(start2), length(end1), length(end2), length(strand))
    if (length(unique(arglen)) != 1 ){
        stop('Arguments are not of the same length')
    }
    dists<- vector()
    for(i in 1:length(start1)){
        if (strand[i] %in% c('+', '-') == FALSE){
            stop(paste('Unexpected strand at line: ', i, ' (', strand[i], ')', sep= ''))
        }
        if (end1[i] < start2[i]) {
            # |---1----|   |----2----|
            x<- (start2[i] - end1[i])
            if (strand[i] == '+') {
                dists[i]<- x
            }
            else {
                dists[i]<- -x ## On negative strand, feature 2 > 1 means upstream
            }
        }
        else if (end2[i] < start1[i]) {
            # |---2----|   |----1----|
            x<- (end2[i] - start1[i])
            if (strand[i] == '+') {
                dists[i]<- x
            }
            else {
                dists[i]<- -x
            }
        }
        else if ((start1[i] >= start2[i] & start1[i] <= end2[i]) | (end1[i] >= start2[i] & end1[i] <= end2[i]) | (start2[i] >= start1[i] & start2[i] <= end1[i]) | (end2[i] >= start1[i] & end2[i] <= end1[i])){
            # |----1----|              |-----1-----|
            #       |-----2----|    |-----2-----|
            dists[i]<- 0
        }
        else {
            stop(paste('Unexpected distance at line', i, '[', start1[i], end1[i], start2[i], end2[i], strand[i], ']'))
        }
    }
    return(dists)
}