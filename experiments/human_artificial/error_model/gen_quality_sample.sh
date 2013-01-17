#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <subsample_size>"
    exit
else
    subsize="$1"
fi

for file in *.fastq
do
    if [ ! -f ${file}.qsample ]
    then
        ### get number of lines in fastq
        reads=`cat $file | sed -n 4~4p | wc -l` 
        ### subsample
        echo "Subsampling $subsize reads from $reads total reads in $file"
        cat $file | sed -n 4~4p | python gen_subsample.py $reads $subsize > $file.qsample
    fi
done

all_size=`cat *.qsample | wc -l`
cat *.qsample | python gen_subsample.py $all_size $subsize > quality.sample
