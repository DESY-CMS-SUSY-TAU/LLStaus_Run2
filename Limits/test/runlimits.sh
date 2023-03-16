#!/bin/bash

CARDDIR=$1

if [ -z "$CARDDIR" ]; then
    echo "\$var is empty"
    exit
fi


for f in $(find $CARDDIR | grep llstau_tauhtauh_2018.txt | sort -V); do
    printf "\n\n\n***** $f *****\n"
    combine -M AsymptoticLimits -t -1 $f
done
