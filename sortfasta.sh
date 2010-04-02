#!/bin/bash

awk 'BEGIN {RS = "\n>ENST"} {gsub(/\n/, "~~", $0); print ">ENST" $0}' $1 | sort | sed 's/~~/\n/g' > $1.sorted;

