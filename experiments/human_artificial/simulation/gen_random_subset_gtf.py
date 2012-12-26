""" This script gets a gffe file and selects a random subset of genes from the specified chromosomes"""
import sys
import random
import pdb

def check_target_lines(fn_in, chrms, coding_only):
    """ Checks the number of valid target lines"""

    target_list = []
    for line in open(fn_in, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        if (len(chrms) > 0) and (not sl[0] in chrms):
            continue
        if coding_only and sl[1] != "protein_coding":
            continue
        tags = sl[8].split('"')
        target_list.append(tags[1])

    return list(set(target_list))

    

def main():
    """Main procedure ..."""

    if len(sys.argv) < 3:
        print >> sys.stderr, "usage: %s <infile> <lines_to_take> [<chr1[,chr2,...]>] ['coding_only']" % sys.argv[0]
        sys.exit(-1)
    
    coding_only = False
    if len(sys.argv) > 4:
        if sys.argv[4] == "coding_only":
            coding_only = True

    if len(sys.argv) > 3:
        chrms = sys.argv[3].split(",")
    else:
        chrms = []

    take_lines = int(sys.argv[2])
    fn_in = sys.argv[1]

    target_list = check_target_lines(fn_in, chrms, coding_only)

    if len(target_list) == 0:
        print >> sys.stderr, "No gene entry for target chromosomes found! - Bailing out!"
        sys.exit(-1)
    elif len(target_list) < take_lines:
        print >> sys.stderr, "Warning - found only %i target lines - take all" % len(target_list)
    else:
        print >> sys.stderr, "Found %i target lines" % len(target_list)

    random.seed(17)
    random.shuffle(target_list)
    target_list = target_list[:min(len(target_list), take_lines)]

    header = True
    took_previous = False
    last_tag = ''
    for line in open(fn_in, 'r'):

        if line[0] == '#':
            if header:
                print line,
            continue

        header = False
        sl = line.strip().split('\t')

        if (len(chrms) > 0 ) and (not sl[0] in chrms):
            continue

        if coding_only and sl[1] != "protein_coding":
            continue

        tags = sl[8].split('"')
        if (tags[1] == last_tag and took_previous) or (tags[1] in target_list):
            print line,
            took_previous = True
        else:
            took_previous = False
        last_tag = tags[1]


    ### check number of target lines
if __name__ == "__main__":
    main()
