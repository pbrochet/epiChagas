#!/bin/bash

INPUT=$1
OUTPUTDIR=$2
REMAP=$3


while read p; do
    echo $p
    echo $REMAP
    grep -w $p $REMAP | awk -F"\t" '{OFS="\t"; print $1, $2, $3}' > $OUTPUTDIR/$p\.bed;
done < $INPUT

