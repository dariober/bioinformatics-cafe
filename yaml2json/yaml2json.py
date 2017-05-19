#!/usr/bin/env python3

docstring= """DESCRIPTION
Convert yaml to json. Read from stdin and write to stdout.
Note that some yaml data types cannot be converted to json.
From 
https://gist.github.com/noahcoad/51934724e0896184a2340217b383af73
USAGE:
cat in.yaml | yaml2.json.py > out.json

Version 0.1.0
"""

import yaml, json, sys

if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help']:
    sys.stderr.write(docstring)
    sys.exit(1)

if len(sys.argv) == 1 or sys.argv[1] == '-':
    fin= sys.stdin
else:
    fin= open(sys.argv[1])

sys.stdout.write(json.dumps(yaml.load(fin), sort_keys=True, indent=2) + '\n')

