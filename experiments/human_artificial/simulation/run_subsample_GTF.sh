#!/bin/bash

set -e 

if [ -z "$1" ]
then
    echo "usage: $0 <lines to subsample> <chr[,chr2,...]> ['coding_only']"
    exit 2
else
    lines=$1
fi
shift

if [ -z "$1" ]
then
    echo "usage: $0 <lines to subsample> <chr[,chr2,...]> ['coding_only']"
    exit 2
else
    chrms="$1"
fi
shift
coding="$1"

work_dir="/fml/ag-raetsch/nobackup2/projects/mmr/human_simulation/annotation"

python gen_random_subset_gtf.py ${work_dir}/hg19.gtf $chrms $lines $coding > ${work_dir}/hg19_`echo $chrms | tr ',' '_'`_subsample_${lines}_genes.gtf 
