#!/bin/bash

set -e

usage() {

   echo "Usage: <genes> <readnum> [<noise>]"
   exit 1
}

[[ -z "$1" ]] && usage
genes=$1
shift

[[ -z "$1" ]] && usage
readnum=$1
shift

noise=""
[[ ! -z "$1" ]] && noise=.noise$1

outdir=${genes}_genes_${readnum}_reads/bowtie/hg19_subsample_${genes}_genes$noise

if [ ! -f ${outdir}/alignment.counts ]
then
    samtools view ${outdir}/alignment.bam | cut -f 3 | sort  -T `pwd`/tmp | uniq -c | sed -e "s/^ *//g" | cut -f 1,2 -d ' ' > ${outdir}/alignment.counts
fi

if [ ! -f ${outdir}/alignment.best.counts ]
then
    samtools view ${outdir}/alignment.bam -F 256 | cut -f 3 | sort -T `pwd`/tmp | uniq -c | sed -e "s/^ *//g" | cut -f 1,2 -d ' ' > ${outdir}/alignment.best.counts
fi
