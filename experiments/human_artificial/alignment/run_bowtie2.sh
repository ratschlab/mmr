#!/bin/bash

set -e

usage() {
    echo "Usage: $0 <genes> <readnum> [<noise>]"
    exit 1
}

[[ -z "$1" ]] && usage
genes=$1
shift

[[ -z "$1" ]] && usage
readnum=$1
shift

noise=""
[[ ! -z "$1" ]] && noise=".noise$1"

chrms=""

bowtie2=/cbio/grlab/share/software/bowtie/bowtie2-2.0.2/bowtie2

transcriptome_idx=/cbio/grlab/nobackup2/projects/mmr/human_simulation_76/annotation/hg19_transcriptome.fa
fastq=${genes}_genes_${readnum}_reads/hg19_${chrms}subsample_${genes}_genes.gtf${noise}.fastq.gz

# generate left and right reads
fastq_left=${fastq%fastq.gz}left.fastq
fastq_right=${fastq%fastq.gz}right.fastq
if [ ! -f $fastq_left ]
then    
    zcat $fastq | grep -e "^@chr.*/1" -A 3 | grep -v -e "^--$" > $fastq_left
fi
if [ ! -f $fastq_right ]
then    
    zcat $fastq | grep -e "^@chr.*/2" -A 3 | grep -v -e "^--$" > $fastq_right
fi

workdir=`pwd`/${genes}_genes_${readnum}_reads
outdir=${workdir}/bowtie/hg19_subsample_${genes}_genes${noise}

mkdir -p $outdir

#$bowtie2 -k 1000 -q -x $transcriptome_idx -1 $fastq_left -2 $fastq_right -S $outdir/alignment.sam 2> $outdir/alignment.log
$bowtie2 -k 200 -q -x $transcriptome_idx -1 $fastq_left -2 $fastq_right -S $outdir/alignment.sam 2> $outdir/alignment.log

samtools view -bS -o ${outdir}/alignment.bam ${outdir}/alignment.sam && rm ${outdir}/alignment.sam

