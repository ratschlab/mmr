#!/bin/bash

set -e

usage() {
    echo "Usage: $0 <genes> <readnum> [noise]"
    exit 1
}

[[ -z "$1" ]] && usage
genes=$1
shift

[[ -z "$1" ]] && usage
reads=$1
shift

noise="$1" 
[[ ! -z "$noise" ]] && noise="noise${noise}."

fname=${genes}_genes_${reads}_reads/tophat/hg19_subsample_${genes}_genes.gtf.${noise}fastq.gz/accepted_hits.bam

[[ ! -f ${fname}.bai ]] && samtools index $fname
[[ ! -f ${fname%bam}best.bam ]] && samtools view -h $fname | python extract_best_tophat.py | samtools view -bS -o ${fname%bam}best.bam -
[[ ! -f ${fname%bam}best.bam.bai ]] && samtools index ${fname%bam}best.bam
