inblast= open('M:/Tritume/best_hits_description.txt', 'r')
outblast= open('M:/Tritume/best_hits_description_brief.txt', 'w')

inblast.readline()
liner= 4
while 1<2:
##    print(liner)
    if liner == 4:
        descr= inblast.readline()
        outblast.write(descr)
        liner= 1
    else:
        descr= inblast.readline()
##        print(descr)
        liner += 1
    if len(descr) == 0:
        break
    
inblast.close()
outblast.close()
    