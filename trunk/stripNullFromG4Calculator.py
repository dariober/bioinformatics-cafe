#!/usr/bin/env python

import sys

docstring= """
DESCRIPTION
    Removes the NULL character '\\x00' from the input file.
    Use this script to remove NULLs from the output of G4calculator.
    
    It appears that G4calculator pads rows with NULLs which cause most programs
    to behave unexpectedly (E.g. Unix "less" will recognize these files as binaries.)

USAGE
    stripNullFromG4Calculator.py g4c.output.txt > g4c.output.fixed.txt
"""

if len(sys.argv) != 2 or sys.argv[1] in ['-h', '--help']:
    sys.exit(docstring)
    
fin= open(sys.argv[1])
for line in fin:
    line= line.strip()
    line= line.replace('\x00', '')
    print(line)
sys.exit()
