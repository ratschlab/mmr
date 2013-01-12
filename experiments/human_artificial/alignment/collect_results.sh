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

work_dir=${genes}_genes_${size}_reads
target=${work_dir}/hg19_${chrms}subsample_${genes}_genes.gtf.${noise}fastq.gz.mapped.${stage}

[[ ! -f $target ]] && samtools merge -n $target ${work_dir}/hg19_${chrms}subsample_${genes}_genes.gtf.${noise}fastq.gz.splits/split_1m.*.gz.mapped.${stage} 

[[ ! -f ${target}.sorted.bam ]] && samtools sort $target ${target}.sorted

[[ ! -f ${target}.sorted.bam.bai ]] && samtools index ${target}.sorted.bam

[[ ! -f ${target}.best.sorted.bam ]] && samtools view -h ${target}.sorted.bam | grep -e "^@" -e "HI:i:0" | samtools view -bS -o ${target}.best.sorted.bam - 

[[ ! -f ${target}.best.sorted.bam.bai ]] && samtools index ${target}.best.sorted.bam

