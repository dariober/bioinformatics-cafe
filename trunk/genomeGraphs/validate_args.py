"""Functions to validate and parse arguments in coverage_screenshots.py
"""
import os

def validate_ymax(ymax):
    """Validate args.ymax
    """
    if ymax == ['max']:
        return(True)
    for y in ymax:
        if y == 'indiv':
            continue
        else:
            try:
                float(y)
            except:
                print('Invalid --ymax: %s' %(y))
                return(False)
    return(True)

def validate_ymin(ymin):
    """Validate args.ymin
    """
    if ymin == ['min']:
        return(True)
    for y in ymin:
        if y == 'indiv':
            continue
        else:
            try:
                float(y)
            except:
                print('Invalid --ymin: %s' %(y))
                return(False)
    return(True)

def parse_names(names, inputlist):
    """Parse list of sample names.
    inputlist:
        Input files (args.ibam) to use when names is None
    """
    if not names:
        snames= [os.path.split(x)[1] for x in inputlist]
        return(snames)
    else:
        return(names)
    