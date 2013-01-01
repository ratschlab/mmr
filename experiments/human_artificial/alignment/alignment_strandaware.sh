#!/bin/bash

echo "Starting alignment_strandaware.sh on machine `uname -a` with args $@"

if [ `hostname` == "node435" ];
then
	echo avoiding `hostname`
	exit 99
fi

export VMEM=`echo \`ulimit -v\` \* 999 / 1000 |bc`
echo consider VMEM=$VMEM
if [[ $VMEM -ne 0 ]];
then
	ulimit -v $VMEM
fi

set -e

export orig_dir=`pwd`

umask 0007
export TERM=linux
export COLORTERM=linux

cd $1 

#echo $1
#echo $2
#echo $3
#echo $4

export stage=$3

if [ $stage == "1" ];
then
	spliced_alignment=0
else
	spliced_alignment=1
fi

#just the left side is considered for the ath alignments
if [ -z "$4" ]
then
	export stranded=0
	echo detected as unstranded 
elif [ $4 == "left" ]
then
	export stranded=1
	export strand="left"
	echo detected as stranded \(left\)
elif [ $4 == "right" ]
then
	export strand="right"
	export stranded=1
	echo detected as stranded \(right\)
fi

source ../../get_tmp_files.sh 

export GENOME=/tmp/rgasp/genomes/hg19_14/hg19.fa #`cat genome.tmp` 
if [ $spliced_alignment == "1" ];
then
	export ACC=/tmp/rgasp/annotations/hg19/acc_pred.bspf #acc.pred.bspf 
	export DON=/tmp/rgasp/annotations/hg19/don_pred.bspf #don.pred.bspf 
fi

umask 0007

export random_wait=`echo "$RANDOM/1000" | bc /dev/stdin`
cd $orig_dir
cd $1 

# determine parameters
if [ $spliced_alignment == "1" ];
then
	export QPALMAPAR=../../parameters.qpalma

	case $stage in 
		2)
            junctionsfile="../../anno.align.8.junctions"
            echo use junctionsfile: $junctionsfile
            export protocol="second"
			export GMPAR=" -M 10 -G 2 -E 10 -l 10 -L 20 -K 12 -C 30 -I 200000 -NI 2 -SA 100 -CT 50 -a -S -seed-hit-cancel-threshold 1000 -report-map-read -report-spliced-read -report-map-region -report-splice-sites 0.95 -filter-max-mismatches 0 -filter-max-gaps 0 -filter-splice-region 5 -qpalma-use-map-max-len 2000 -f bamn -samtools /cbio/grlab/share/software/samtools/ -polytrim 40 -qpalma-prb-offset-fix -min-spliced-segment-len 8 -junction-remapping-coverage 2 -junction-remapping $junctionsfile -score-annotated-splice-sites $junctionsfile -acc $ACC/pred/contig_%i%c -don $DON/pred/contig_%i%c -report-splice-sites-top-perc 0.01 -QMM 7 -max-dp-deletions 3"
			;;
		3)
            junctionsfile="../../anno.align.8.junctions"
            echo use junctionsfile: $junctionsfile
            export protocol="second"
			export GMPAR=" -M 6 -G 2 -E 6 -l 10 -L 20 -K 12 -C 30 -I 200000 -NI 2 -SA 100 -CT 50 -a -S -seed-hit-cancel-threshold 1000 -report-map-read -report-spliced-read -report-map-region -report-splice-sites 0.95 -filter-max-mismatches 0 -filter-max-gaps 0 -filter-splice-region 5 -qpalma-use-map-max-len 2000 -f bamn -samtools /cbio/grlab/share/software/samtools/ -polytrim 40 -qpalma-prb-offset-fix -min-spliced-segment-len 8 -junction-remapping-coverage 2 -junction-remapping $junctionsfile -score-annotated-splice-sites $junctionsfile -acc $ACC/pred/contig_%i%c -don $DON/pred/contig_%i%c -report-splice-sites-top-perc 0.01 -QMM 7 -max-dp-deletions 3"
			;;
		
	esac
	if [ ! -f $QPALMAPAR ]
		then
		echo qpalma parameter file $QPALMAPAR does not exist
		exit -1
	fi
	echo using $QPALMAPAR as qpalma parameters
else
	export GMPAR=" -M 6 -G 2 -E 6 -l 18 -z 10 -seed-hit-cancel-threshold 10000  -report-map-read -report-spliced-read  -report-map-region -report-splice-sites 0.8 -f bamp -rtrim 40 -rtrim-step 6"
fi

if [ -f genomemapper.par ]
then
	echo using the parameters given in $1/genomemapper.par
	export GMPAR=`cat genomemapper.par`
fi
echo using the following parameters for genomemapper: $GMPAR

if [ -z "$NSLOTS" ];
then
	NSLOTS=1 ;
	echo Warning: NSLOTS not set, setting to $NSLOTS
fi
NSLOTS=1

mem=`cat /proc/meminfo  | grep MemTotal | cut -f 2 -d :| cut -f 1 -d k`

if [ $mem -lt 40000000 ]
then
	NSLOTS=1
fi
echo NSLOTS=$NSLOTS

palmapper_path=/cbio/grlab/home/akahles/git/software/palmapper

if [ $spliced_alignment == "1" ];
then
	if [ $stranded == "1" ];
	then
		echo $palmapper_path/palmapper -i $GENOME -q ../../$2 -o ../../$2.mapped.$stage -stranded $strand -protocol $protocol -threads $NSLOTS -qpalma $QPALMAPAR $REPORT $REPORT_RO $GMPAR -index-precache -samtools /cbio/grlab/share/software/samtools/
		#valgrind --max-stackframe=4097344 
		#gdb --args 
		time $palmapper_path/palmapper -i $GENOME -q ../../$2 -o ../../$2.mapped.$stage -stranded $strand -protocol $protocol -threads $NSLOTS -qpalma $QPALMAPAR $REPORT $REPORT_RO $GMPAR  -index-precache -samtools /cbio/grlab/share/software/samtools/
	else 
		echo $palmapper_path/palmapper -i $GENOME -q ../../$2 -o ../../$2.mapped.$stage -threads $NSLOTS -qpalma $QPALMAPAR $REPORT $REPORT_RO $GMPAR  -index-precache -samtools /cbio/grlab/share/software/samtools/
		time $palmapper_path/palmapper -i $GENOME -q ../../$2 -o ../../$2.mapped.$stage -threads $NSLOTS -qpalma $QPALMAPAR $REPORT $REPORT_RO $GMPAR -index-precache -samtools /cbio/grlab/share/software/samtools/
	fi 
else
	if [ $stranded == "1" ];
	then
		echo $palmapper_path/palmapper -i $GENOME -q ../../$2 -o ../../$2.mapped.$stage -threads $NSLOTS $REPORT $REPORT_RO $GMPAR -stranded $strand -protocol $protocol -index-precache -samtools /cbio/grlab/share/software/samtools/
		time $palmapper_path/palmapper -i $GENOME -q ../../$2 -o ../../$2.mapped.$stage -threads $NSLOTS $REPORT $REPORT_RO $GMPAR -stranded $strand -protocol $protocol -index-precache -samtools /cbio/grlab/share/software/samtools/
	else
		echo $palmapper_path/palmapper -i $GENOME -q ../../$2 -o ../../$2.mapped.$stage -threads $NSLOTS $REPORT $REPORT_RO $GMPAR -index-precache -samtools /cbio/grlab/share/software/samtools/
		time $palmapper_path/palmapper -i $GENOME -q ../../$2 -o ../../$2.mapped.$stage -threads $NSLOTS $REPORT $REPORT_RO $GMPAR -index-precache -samtools /cbio/grlab/share/software/samtools/
	fi 
fi

touch ../../$2.done.$stage
