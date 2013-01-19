#!/bin/bash

set -e 

if [ -z "$1" ]
then
    echo "Usage: $0 <BAM-file>"
    exit 1
else
    filename=$1
fi

samtools view -h $filename | python `dirname $0`/extract_multimapper.py | samtools view -bS -o ${filename}.multimappers.bam - 
#[[ "$filename" == "accepted_hits.ID_sorted.bam" ]] && ln -s accepted_hits.ID_sorted.multimappers.bam accepted_hits.multimappers.bam && filename=accepted_hits.multimappers.bam
samtools sort ${filename}.multimappers.bam ${filename}.multimappers.sorted
samtools index ${filename}.multimappers.sorted.bam
