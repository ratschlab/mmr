""" This script gets a gffe file and selects a random subset of genes from the specified chromosomes"""
import sys
import random

def check_target_lines(fn_in, chrms, coding_only):
    """ Checks the number of valid target lines"""

    target_list = []
    line_idx = 0
    for line in open(fn_in, 'r'):
        line_idx += 1
        if line[0] == '#':
            continue
        sl = line.strip().split()
        if not sl[0] in chrms:
            continue
        if not sl[2] == "gene":
            continue
        if coding_only and sl[1] != "protein_coding":
            continue
        target_list.append(line_idx)

    return target_list

    

def main():
    """Main procedure ..."""

    if len(sys.argv) < 4:
        print >> sys.stderr, "usage: %s <infile> <chr[,chr2,...]> <lines_to_take> ['coding_only']" % sys.argv[0]
        sys.exit(-1)
    
    coding_only = False
    if len(sys.argv) > 4:
        if sys.argv[4] == "coding_only":
            coding_only = True

    fn_in = sys.argv[1]
    chrms = sys.argv[2].split(",")
    take_lines = int(sys.argv[3])

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

    line_idx = 0
    we_take_it = False
    for line in open(fn_in, 'r'):
        line_idx += 1

        if line[0] == '#':
            print line,
            continue
        sl = line.strip().split()

        if not sl[0] in chrms:
            continue

        if coding_only and sl[1] != "protein_coding":
            continue

        if sl[2] == "gene":
            we_take_it = (line_idx in target_list)

        if we_take_it:
            print line,
 


    ### check number of target lines
if __name__ == "__main__":
    main()
