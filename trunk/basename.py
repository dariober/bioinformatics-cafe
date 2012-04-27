#!/home/berald01/.local/bin/python

import sys
import re
import os

if len(sys.argv) == 1:
    sys.exit("""
Strip file path and regex from string. Similar to bash basename but using python re.sub() instead

USAGE
    basename.py <string> <regex>

EXAMPLE
    ## Return 'myfile'
    basename.py '\.bam$' path/to/myfile.bam 
""")

x= os.path.split(sys.argv[1])[1]
print(re.sub(sys.argv[2], '', x))
sys.exit()