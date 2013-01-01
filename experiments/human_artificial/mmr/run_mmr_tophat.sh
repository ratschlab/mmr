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
    noise="noise${noise}."
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

mmr="$HOME/git/software/RNAgeeq/mm_resolve/threaded_oop_mip/mmr"
samtools=/cbio/grlab/share/software/samtools/samtools

filename="${genes}_genes_${size}_reads/tophat/hg19_${chrms}subsample_${genes}_genes.gtf.${noise}fastq.gz/accepted_hits.ID_sorted.bam"
target=${filename}.mmr${filter}
if [ ! -f "$target" ]
then
    time $mmr -w 20 -t $threads -p -b -T $samtools -f -F $filter -v -I 3 -o $target $filename > ${target}.log
    $samtools sort $target ${target}.sorted
    $samtools index ${target}.sorted.bam
else
    echo "$target already exists"
fi
