#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "usage: $0 <genes> <readnum> [noise] [<chr1[,chr2,...]>]"
    exit
else
    genes=$1
fi
shift

if [ -z "$1" ]
then
    echo "usage: $0 <genes> <readnum> [noise] [<chr1[,chr2,...]>]"
    exit
else
    size=$1
fi

noise="$2"
if [ ! -z "$noise" ]
then
    noise=noise${noise}. 
fi

chrms="$3"
if [ ! -z "$chrms" ]
then
    chrms=`echo $chrms | tr ',' '_'`_
fi

fastq=${genes}_genes_${size}_reads/hg19_${chrms}subsample_${genes}_genes.gtf.${noise}fastq.gz

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

### set paths
export PATH=/cbio/grlab/share/software/bowtie-0.12.7:$PATH
export PATH=/cbio/grlab/share/software/bowtie2-2.0.2:$PATH
export PATH=/cbio/grlab/share/software/samtools:$PATH

# align reads
outdir=${genes}_genes_${size}_reads/tophat/hg19_${chrms}subsample_${genes}_genes.gtf.${noise}fastq.gz
anno=/cbio/grlab/nobackup/projects/rgasp/annotations/human/hg19/encode_data.rel2b.liftover.all.gtf

### parameters
mismatches=6
gaps=2
edit=6
threads=8
max_hits=50

mkdir -p $outdir
tophat2=/cbio/grlab/share/software/tophat-2.0.5/bin/tophat2

$tophat2 /cbio/grlab/nobackup/projects/rgasp/genomes/hg19_14/hg19.fa $fastq_left $fastq_right --output-dir $outdir --read-mismatches $mismatches --read-edit-dist $edit --read-gap-length $gaps --max-intron-length 50000 --max-multhits $max_hits --num-threads $threads --GTF $anno --report-secondary-alignments 2> $outdir/run_tophat.log
