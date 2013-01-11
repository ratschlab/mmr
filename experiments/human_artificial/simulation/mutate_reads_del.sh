#!/bin/bash

set -e 

if [ -z "$1" ]
then
    echo "usage: <genes> <readnum> [<noise>] [<chr1,[chr2,...]>]"
    exit 2
else
    genes="$1"
fi
shift

if [ -z "$1" ]
then
    echo "usage: <genes> <readnum> [<noise>] [<chr1,[chr2,...]>]"
    exit 2
else
    readnum="$1"
fi
noise="$2"
noisetag=""
if [ ! -z "$noise" ]
then
    noisetag=".noise${noise}"
    noise="-n $noise" 
fi

chrms="$3"
if [ ! -z "$chrms" ]
then
    chrms="`echo $chrms | tr ',' '_'`_"
fi

#sample="$HOME/git/projects/2011/rgasp/GenArtificialReads/sample_from_flux.py"
sample="$HOME/git/projects/2011/rgasp/GenArtificialReads/sample.py"
workdir=`pwd`

export PYTHONPATH=$HOME/git/tools/python:$PYTHONPATH
filebase=${genes}_genes_${readnum}_reads/hg19_${chrms}subsample_${genes}_genes.gtf
gio_file=/cbio/grlab/nobackup/projects/rgasp/genomes/hg19_14/hg19.gio/genome.config

python $sample -s -F -q ${workdir}/error_model/quality.sample $noise -G $gio_file -m ${workdir}/error_model/mismatch.pickle -o ${filebase}${noisetag}.optimal -f ${filebase}.fasta -b ${filebase}.bed | gzip -c -6 > "${filebase}${noisetag}".fastq.gz
