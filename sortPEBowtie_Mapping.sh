#!/bin/bash

awk ' BEGIN { cntr = 1 }
{
 if ( cntr == 1 ) {
	printf("%s~~", $0)
 } 
 if ( cntr == 2 ) {
	printf("%s", $0)
 }
  if ( cntr == 2 ) {
    printf("\n")
    cntr = 0
  }
  cntr++
}
END {
  rest = int(cntr)
  if ( rest != 1 ) {
    printf("\n")
  }
}
' $1 | sort -T . -S 4G -k 2,3 | sed -e 's/~~/\n/g' -e 's/\n\n/\n/g' > $1.sorted
