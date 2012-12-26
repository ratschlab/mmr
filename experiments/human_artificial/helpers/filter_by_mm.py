"""This script takes an alignment file ins SAM format from
stdin and filter according to the given max # of edit ops."""
import sys

if len(sys.argv) < 2:
    print >> sys.stderr, "Usage: %s <max_mm>" % sys.argv[0]
    sys.exit(1)

max_mm = int(sys.argv[1])

for line in sys.stdin:
    if line[0] == '@':
        print line,
        continue

    sl = line.strip().split()
    for i in sl[11:]:
        if i[:2] == 'NM':
            if int(i[5:]) <= max_mm:
                print line,
