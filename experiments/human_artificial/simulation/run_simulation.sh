#!/bin/bash
# (C)   Andre Kahles

set -e

usage() {
    echo "usage: $0 <genes> <readnum> <readlen> [<chr1[,chr2,...]>]"
    exit
}

[[ -z "$1" ]] && usage
genes="$1"
shift

[[ -z "$1" ]] && usage
size="$1"
shift

[[ -z "$1" ]] && usage
readlen="$1"
shift

chrms="$1"

flux=/cbio/grlab/share/software/FluxSimulator/flux-simulator-1.1.1-20121103021450/bin/flux-simulator

workdir="/cbio/grlab/nobackup2/projects/mmr/human_simulation_${readlen}"
if [ ! -z "$chrms" ]
then
    gtf=$workdir/annotation/hg19_`echo $chrms | tr ',' '_'`_subsample_${genes}_genes.gtf
else
    gtf=$workdir/annotation/hg19_subsample_${genes}_genes.gtf
fi
flux_par_file=${workdir}/flux_parameter_${genes}_genes_${size}_reads.par
outdir=${workdir}/${genes}_genes_${size}_reads
mkdir -p $outdir

cp flux_parameter_stub $flux_par_file

MOL_NUM=$((4*$size))
READ_NUM=$size
echo "NB_MOLECULES    $MOL_NUM" >> $flux_par_file
echo "READ_NUMBER   $READ_NUM" >> $flux_par_file
echo "READ_LENGTH   $readlen" >> $flux_par_file
echo "REF_FILE_NAME $gtf" >> $flux_par_file
#if [ ! -f ${workdir}/76_error.model ]
#then
#    echo "$workdir does not contain the file 76_error.model" 
#    exit -1
#fi
#echo "ERR_FILE      ${workdir}/76_error.model" >> $flux_par_file

gtf_base=`basename $gtf`
echo "LIB_FILE_NAME ${outdir}/${gtf_base}.lib" >> $flux_par_file
echo "PRO_FILE_NAME ${outdir}/${gtf_base}.pro" >> $flux_par_file
echo "SEQ_FILE_NAME ${outdir}/${gtf_base}.bed" >> $flux_par_file

export TMPDIR=${workdir}/tmp

mkdir -p $TMPDIR
mkdir -p ${workdir}/log_flux

$flux --log DEBUG --force -p $flux_par_file &> ${workdir}/log_flux/flux_sim_${genes}_genes_${size}_reads.log

cat ${outdir}/${gtf_base}.bed | grep -v -e polyA > ${outdir}/${gtf_base}.noPolyA.bed

bedToBam=/cbio/grlab/share/software/BEDTools/BEDTools-Version-GIT/bin/bedToBam
genome=${workdir}/annotation/hg19.genome

$bedToBam -bed12 -i ${outdir}/${gtf_base}.bed -g $genome > ${outdir}/${gtf_base}.bam
$bedToBam -bed12 -i ${outdir}/${gtf_base}.noPolyA.bed -g $genome > ${outdir}/${gtf_base}.noPolyA.bam
echo done
