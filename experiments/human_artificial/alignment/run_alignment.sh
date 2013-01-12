#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "usage: $0 <genes> <readnum> <stage> [noise] [<chr1[,chr2,...]>]"
    exit
else
    genes=$1
fi
shift

if [ -z "$1" ]
then
    echo "usage: $0 <genes> <readnum> <stage> [noise] [<chr1[,chr2,...]>]"
    exit
else
    size=$1
fi
shift

if [ -z "$1" ]
then
    echo "usage: $0 <genes> <readnum> <stage> [noise] [<chr1[,chr2,...]>]"
    exit
else
    stage=$1
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


if [ ! -f ${fastq}.splits/split_1m.000 -a ! -f ${fastq}.splits/split_1m.000.gz ];
then
    mkdir -p ${fastq}.splits
    zcat $fastq | split --suffix-length=3 -d -l 1000000 - ${fastq}.splits/split_1m.
    for file in `find ${fastq}.splits/split_1m.[0-9][0-9][0-9]`
    do
        gzip -6 $file
    done
fi

for file in `find ${fastq}.splits/split_1m.[0-9][0-9][0-9].gz`
do
    if [ ! -f ${file}.done.$stage ]
    then
        echo submitting $file
        echo time ./alignment_strandaware.sh ${fastq}.splits/ $file $stage $strand | qsub -o $file.log.$stage -j y -l h_vmem=30G -cwd -N hu.$stage -pe make 2 >& $file.jobid.$stage
    fi
done


