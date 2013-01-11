#!/bin/bash

set -e 

if [ -z "$1" ]
then
    echo "usage: $0 <genes to subsample> <chr[,chr2,...]> ['coding_only']"
    exit 2
else
    genes=$1
fi

chrms="$2"
coding="$3"

work_dir="/cbio/grlab/nobackup2/projects/mmr/human_simulation/annotation"

if [ ! -z "$chrms" ]
then
    python gen_random_subset_gtf.py ${work_dir}/hg19.gtf $genes $chrms $coding > ${work_dir}/hg19_`echo $chrms | tr ',' '_'`_subsample_${genes}_genes.gtf 
else
    python gen_random_subset_gtf.py ${work_dir}/hg19.gtf $genes $chrms $coding > ${work_dir}/hg19_subsample_${genes}_genes.gtf 
fi
