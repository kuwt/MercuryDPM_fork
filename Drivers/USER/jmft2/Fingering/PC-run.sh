#!/bin/bash
#
# Run the simulation and the post-processing script.

set -e

if [[ $# == 2 && $2 == '-r' ]]; then
    nice -n 19 ./PerryComb $1/Exp.config -r >> $1/Exp.txt 2>> $1/Exp.cerr
elif [ $# == 1 ]; then
    nice -n 19 ./PerryComb $1/Exp.config > $1/Exp.txt 2> $1/Exp.cerr
else
    printf "Usage: $0 Exp [-r]\n"
    exit 1
fi

./PC-post.sh $1
