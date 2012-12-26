import sys

print >> sys.stdout, 'gene\ttranscript\tFPKM\tcov'
for line in open(sys.argv[1], 'r'):
    sl = line.strip().split('\t')
    if sl[2] != 'transcript':
        continue
    tags = sl[8].split(';')
    gene_id = tags[0].split()
    tra_id = tags[1].split()
    fpkm = tags[2].split()
    cov = tags[6].split()
    print >> sys.stdout, '%s\t%s\t%s\t%s' %(gene_id[1][1:-1], tra_id[1][1:-1], fpkm[1][1:-1], cov[1][1:-1])
