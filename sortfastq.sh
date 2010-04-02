#!/bin/bash

awk 'BEGIN {RS = "\n@ENST"} {gsub(/\n/, "~~", $0); print "@ENST" $0}' $1 | sort -S 4G | sed -e 's/~~/\n/g' -e 's/\n\n/\n/g'  > $1.sorted;
