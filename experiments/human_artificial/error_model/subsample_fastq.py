#!/usr/bin/python

# written by Vipin T Sreedharan and Gunnar Raetsch

# Usage : 
# 1. Provide a valid fastq file with reads.
# 2. Number of reads to be extracted from the file.
# 3. Provide result file in fastq format.

import sys
from Bio import SeqIO
import random 

# Accept a fastq file from user
in_file = sys.argv[1]
# Number of reads to be extracted from the fastq file 
target_no = int(sys.argv[2])
# Result file 
out_file = sys.argv[3]

print("")
print("subsample_fastq.py version 0.1")
print("------------------------------")
#print("  counting fastq items ... "),

# open the file 
in_handle = open(in_file, "rU")
no_reads = 0 
for rec in SeqIO.parse(in_handle, "fastq"):
    no_reads += 1
in_handle.close()
#print("  done") 

try:
  accept_prob = (1.0*target_no)/no_reads
except:
  print("\n  Computing acceptance probability failed, setting to 1\n") 
  accept_prob = 1 

# Input verification 
if target_no > no_reads:
    print("\n  Warning:\tLess number of reads present in fastq file than target number of reads to extract.")
    print("\t\tTotal number of reads in fastq file : " + str(no_reads))
    print("\t\tNumber of reads to be extracted : " + str(target_no) + "\n")
    accept_prob = 1 
	
# Taking random reads.    
#print ''
#print "-----Summary of random fastq filtering-----"
print "  Total number of short reads present in the input file : ", no_reads
print "  Probability of accepting a read:                        ", accept_prob
in_handle = open(in_file, "rU")
res_read = 0
c = 0 
res_file = open(out_file, "w+")
for rec in SeqIO.parse(in_handle, "fastq"):
    random_number = random.random()
    c += 1 
    if random_number <=  accept_prob:
        res_read += 1
        res_file.write(rec.format("fastq"))
    if target_no == res_read: # stop when we have enough reads
        break
print "  Number of reads taken :                                 ", res_read
print "  Total number of reads scanned :                         ", c
#print '----------------------------------------------'
res_file.close()
in_handle.close()
