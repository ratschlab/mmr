#!/bin/bash

set -e 

usage() {
    echo "usage: $0 <fastq-sample> <ref-seq> <quality_sample>"
    exit 2
}

[[ -z "$1" ]] && usage
fastq="$1"
shift
[[ -z "$1" ]] && usage
ref_seq="$1"
shift
[[ -z "$1" ]] && usage
quality="$1"

error_model=mismatch.pickle

path_to_palmapper="$HOME/git/software/palmapper/"

aligned_reads=`mktemp`

### create indices
if [ ! -f $ref_seq.cid ]
then
	echo ""
    echo "creating index for $ref_seq ..."
    cp $ref_seq $gmapper_index_dir/
    $path_to_palmapper/pmindex -v -i $gmapper_index_dir/`basename $ref_seq`
    random_index_dir=1
else
    gmapper_index_dir=`dirname $ref_seq`
fi

### align
echo "starting alignment"
echo ""

echo $path_to_palmapper/palmapper -q $fastq -i $gmapper_index_dir/`basename $ref_seq` -o $aligned_reads -M 6 -G 2 -E 6 -f bedx -z 5 -seed-hit-cancel-threshold 10000 -index-precache
$path_to_palmapper/palmapper -q $fastq -i $gmapper_index_dir/`basename $ref_seq` -o $aligned_reads -M 6 -G 2 -E 6 -f bedx -z 5 -seed-hit-cancel-threshold 10000 -index-precache

echo ""
echo "generate mismatch statistics"
echo ""
echo "python $HOME/git/projects/2011/rgasp/GenArtificialReads/tools/mismatch/mismatch.py $aligned_reads $error_model 8"
python $HOME/git/projects/2011/rgasp/GenArtificialReads/tools/mismatch/mismatch.py $aligned_reads $error_model 8

echo "cleaning up"
echo ""
rm $aligned_reads


