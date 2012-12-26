import sys
import re

pat = re.compile(r"HI:i:[1-9]")

for line in sys.stdin:
    if line[0] == '@':
        print line,
        continue
    if pat.search(line):
        continue
    print line,

