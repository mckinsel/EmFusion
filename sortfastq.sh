#!/bin/bash


awk ' BEGIN { cntr = 1 }
{
	if ( cntr < 4 ) {
		printf("%s~~", $0)
	}
	if ( cntr == 4 ) {
		printf("%s\n", $0)
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
' $1 | sort -T . -S 4G | sed -e 's/~~/\n/g' -e 's/\n\n/\n/g'  > $1.sorted;
