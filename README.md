MMR: A Tool for Read Multi-Mapper Resolution 
============================================

Motivation:
-----------
Mapping high throughput sequencing data to a reference genome is an
essential step for most analysis pipelines aiming at the computational
analysis of genome and transcriptome sequencing data. Breaking ties
between equally well mapping locations poses a severe problem not only
during the alignment phase, but also has significant impact on the
results of downstream analyses. We present the multimapper resolution
(MMR) tool that infers optimal mapping locations from the coverage
densities of other mapped reads	in vicinity.

Results:
--------
Filtering alignments with MMR can significantly improve the
performance of downstream analyses like transcript quantitation and
differential testing. We illustrate that the accuracy (Spearman
correlation) of transcript quantification increases by 17% when
using reads of length 51. In addition, MMR decreases the alignment
file sizes and this leads to a reduction of running time of the
quantification tool.  Our efficient implementation of the MMR
algorithm is easily applicable as a post-processing step to existing
alignment files in BAM format. Its complexity scales linearly with the
number of alignments and requires no further inputs.

Availability:
-------------
The source code is public under GPLv3 and can be downloaded from 
GitHub (https://github.com/ratschlab/mmr).
Supplemental text and figures, comprehensive testing results and
further information can be found at http://bioweb.me/mmr.

Documentation:
--------------
The user documentation and a description of the commandline parameters 
are available in file DOCUMENTATION.

