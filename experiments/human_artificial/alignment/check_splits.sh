#!/bin/bash

set -e

usage() {
    echo "Usage: $0 <stage> ['cleanup']"
    exit 1
}

[[ -z "$1" ]] && usage
stage="$1"
shift

what="$1"
counter_uf=0
counter_fi=0


for file in *.log.2;
do
    if [ -z "`cat ${file%log.$stage}log.2 | grep -e \"samtools subprocess terminated successfully\"`" ]
    then 
        echo $file
        counter_uf=$(($counter_uf + 1))
        if [ "$what" == "cleanup" ]
        then
            rm -f $file
            rm -f ${file%log.$stage}done.$stage
        fi  
    else
        counter_fi=$(($counter_fi + 1))
    fi
done
echo unfinished: $counter_uf
echo finished: $counter_fi
