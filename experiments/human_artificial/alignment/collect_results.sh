#!/bin/bash

set -e

usage() {
    echo "usage: $0 <genes> <readnum> <stage> [noise] [<chr1,[chr2,...]>]"
}

[[ -z "$1" ]] && usage && exit 2
genes="$1"
shift
[[ -z "$1" ]] && usage && exit 2
size="$1"
shift
[[ -z "$1" ]] && usage && exit 2
stage="$1"
shift

noise="$1"
if [ ! -z "$noise" ]
then
    noise="noise${noise}."
fi

chrms="$2"
if [ ! -z "$chrms" ]
then
    chrms="`echo $chrms | tr ',' '_'`_"
fi

samtools=/cbio/grlab/share/software/samtools/samtools
work_dir=${genes}_genes_${size}_reads

$samtools merge -n ${work_dir}/hg19_${chrms}subsample_${genes}_genes.gtf.${noise}fastq.gz.mapped.${stage} ${work_dir}/hg19_${chrms}subsample_${genes}_genes.gtf.${noise}fastq.gz.splits/split_1m.*.gz.mapped.${stage} 

