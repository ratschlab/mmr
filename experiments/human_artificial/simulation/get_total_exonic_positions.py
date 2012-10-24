import sys

if len(sys.argv) < 2:
    print >> sys.stderr, "usage: %s <GFF3>" % sys.argv[0]

exon_pos = 0
for line in open(sys.argv[1], 'r'):
    if line[0] == '#':
        continue
    sl = line.strip().split()
    if sl[2] == 'exon':
       exon_pos += (int(sl[4]) - int(sl[3])) 

print >> sys.stdout, "file %s has %i exonic positions" % (sys.argv[1], exon_pos)
