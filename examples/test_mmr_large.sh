#!/bin/bash

set -e

which samtools > /dev/null 2>&1 || (echo -e "samtools has not been found to be present in your PATH but is required to run this test. Please add the location of samtools via:\n\nexport PATH=<path_to_samtools>:\$PATH"; exit 1;)

testfile=example_file_large.ID_sorted.bam

outfile=${testfile%bam}mmr.bam
echo -e "\nRunning: ../mmr -F 0 -p -b -v -o $outfile $testfile \n"
../mmr -F 0 -p -b -v -o $outfile $testfile

mkfifo example_file_large.ID_sorted.mmr.sam
samtools view example_file_large.ID_sorted.mmr.bam > example_file_large.ID_sorted.mmr.sam &
md5sum -c --quiet ${outfile}.md5 && echo -e "\nTest successfully passed" || echo -e "\nTest failed"
rm example_file_large.ID_sorted.mmr.sam
    
