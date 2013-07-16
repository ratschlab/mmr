#!/bin/bash

set -e

usage() {
    echo "usage: $0 <genes> <size> <filter> [<noise>] [<threads>] [chrms]"
    exit 2
}


[[ -z "$1" ]] && usage
genes="$1"
shift

[[ -z "$1" ]] && usage
size="$1"
shift

[[ -z "$1" ]] && usage
filter="$1"
shift

noise=$1
if [ "$noise" == "-" ]
then
    noise=""
fi
if [ ! -z "$noise" ]
then
    echo $noise
    noise=".noise${noise}"
fi

threads="1"
if [ ! -z "$2" ]
then
    threads="$2"
fi
chrms=$3
if [ ! -z "$chrms" ]
then
    chrms=`echo $chrms | tr ',' '_'`_
fi

mmr="$HOME/git/software/RNAgeeq/mm_resolve/mmr"

filename="${genes}_genes_${size}_reads/bowtie/hg19_${chrms}subsample_${genes}_genes${noise}/alignment.ID_sorted.bam"
if [ ! -f "$filename" ]
then    
    echo "sorting ${filename%ID_sorted.bam}bam"
    samtools sort -n ${filename%ID_sorted.bam}bam ${filename%.bam}
fi

target=${filename}.mmr${filter}
if [ ! -f "$target" ]
then
    #echo $mmr -r 2 -w 20 -t $threads -p -b -T samtools -F $filter -v -I 5 -o $target $filename
    time $mmr -r 2 -w 20 -t $threads -p -b -T samtools -F $filter -v -I 3 -A 500 -o $target $filename > ${target}.log
    #time $mmr -r 2 -w 30 -t $threads -p -b -T samtools -F $filter -v -I 3 -o $target $filename > ${target}.log
    #time $mmr -r 2 -w 12 -t $threads -p -b -T samtools -F $filter -v -I 3 -o $target $filename > ${target}.log
else
    echo "$target already exists"
fi

if [ ! -f ${target}.counts ]
then
    samtools view $target | cut -f 3 | sort | uniq -c | sed -e "s/^ *//g" | cut -f 1,2 -d ' ' > ${target}.counts
    oldpwd=$(pwd)
    cd $(dirname $target)
    ln -s $(basename $target).counts $(echo $(basename $target) | sed -e "s/.ID_sorted.bam//g").counts
    cd $oldpwd
else
    echo ${target}.counts already counted
fi

