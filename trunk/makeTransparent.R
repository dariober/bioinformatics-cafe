makeTransparent<-function(someColor, alpha=100){
  "Given a colour name (e.g. 'red'), make it transparent. 
  someColor:
        Vector of colour names to make transparent e.g. c('red', 'blue')
  alpha:
        Alpha transparency. 100 fully opaque, 0 fully transparent.
  Credit: http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
  "
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}