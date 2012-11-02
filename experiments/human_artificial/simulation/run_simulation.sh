#!/bin/bash
# (C)   Andre Kahles

set -e

if [ -z "$1" ]
then
    echo "usage: $0 <gtf> <size>"
    exit
else
    gtf="$1"
fi
shift

if [ -z "$1" ]
then
    echo "usage: $0 <gtf> <size>"
    exit
else
    size="$1"
fi

flux=/fml/ag-raetsch/share/software/flux-simulator-1.1/bin/flux-simulator
workdir="/fml/ag-raetsch/nobackup2/projects/mmr/human_simulation"
flux_par_file=${workdir}/flux_parameter_${size}_reads.par
outdir=${workdir}/${size}_reads
mkdir -p $outdir

cp flux_parameter_stub $flux_par_file

MOL_NUM=$((8*$size))
READ_NUM=$((2*$size))
echo "NB_MOLECULES    $MOL_NUM" >> $flux_par_file
echo "READ_NUMBER   $READ_NUM" >> $flux_par_file
echo "REF_FILE_NAME $gtf" >> $flux_par_file
if [ ! -f ${workdir}/76_error.model ]
then
    echo "$workdir does not contain the file 76_error.model" 
    exit -1
fi
gtf_base=`basename $gtf`
echo "ERR_FILE      ${workdir}/76_error.model" >> $flux_par_file
echo "LIB_FILE_NAME ${outdir}/${gtf_base}.lib" >> $flux_par_file
echo "PRO_FILE_NAME ${outdir}/${gtf_base}.pro" >> $flux_par_file
echo "SEQ_FILE_NAME ${outdir}/${gtf_base}.bed" >> $flux_par_file

$flux --log DEBUG --force -p $flux_par_file &> ${workdir}/flux_sim_${size}_reads.log

bedToBam=/fml/ag-raetsch/share/software/BEDTools-Version-2.16.2/bin/bedToBam
genome=${workdir}/annotation/hg19.genome

$bedToBam -bed12 -i ${outdir}/${gtf_base}.bed -g $genome > ${outdir}/${gtf_base}.bam

cat ${outdir}/${gtf_base}.bed | grep -v -e polyA > ${outdir}/${gtf_base}.noPolyA.bed
$bedToBam -bed12 -i ${outdir}/${gtf_base}.noPolyA.bed -g $genome > ${outdir}/${gtf_base}.noPolyA.bam
echo done
