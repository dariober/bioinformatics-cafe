#!/bin/bash

./yaml2json.py test_data/in.yaml

cat test_data/in.yaml | ./yaml2json.py

cat test_data/in.yaml | ./yaml2json.py -

./yaml2json.py -h
