#!/bin/bash

set -e

testfile=example_file_small.ID_sorted.bam

outfile=${testfile%bam}mmr.bam
echo -e "\nRunning: ../mmr -F 0 -b -v -o $outfile $testfile \n"
../mmr -F 0 -b -v -o $outfile $testfile

md5sum -c --quiet ${outfile}.md5 && echo -e "\nTest successfully passed" || echo -e "\nTest failed"
    
