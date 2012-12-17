"""This script reads an alignment file in SAM format from standard in and returns all alignments of reads, that have more than one mapping location"""

import sys 

last_id = ''

curr_lines = []

for line in sys.stdin:
    
    if line[0] == '@':
        print line,
        continue

    sl = line.strip().split()

    if last_id == '':
        last_id = sl[0]

    if last_id == sl[0]:
        curr_lines.append(line)
    else:
        if len(curr_lines) > 1:
            for ll in curr_lines:
                print ll,
        curr_lines = [line]
        last_id = sl[0]
        

if len(curr_lines) > 1:
    for ll in curr_lines:
        print ll,
    
    
    
