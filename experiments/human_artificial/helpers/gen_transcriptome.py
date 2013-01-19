import sys
import pdb
import genome_utils

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--annotation', dest='gtf', metavar='FILE', help='annotation in GTF format', default='-')
    required.add_option('-g', '--genome', dest='fasta', metavar='FILE', help='genome sequence in FASTA format', default='-')
    required.add_option('-o', '--outfile', dest='outfile', metavar='FILE', help='outfile to store transcriptome in FASTA format', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()

    if len(argv) < 3 or options.gtf == '-' or options.fasta == '-':
        parser.print_help()
        sys.exit(1)

    return options

def main():
    """The main routine of this Script"""

    options = parse_options(sys.argv)
    log = sys.stdout
    trans_dict = dict()
    options.chr_prefix = 'chr'

    ### parse gtf file once to get the transcript coordinates
    print >> log, "Parsing transcripts from %s ..." % options.gtf
    counter = 0
    for line in open(options.gtf, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        if sl[2] != "exon":
            continue

        if options.verbose and counter % 10000 == 0:
            print >> log, "[ %i exons ]" % counter
        counter += 1
        tags = sl[8].split(';')
        trans_id = tags[1].split('"')
        trans_id = trans_id[1]
        try:
            if (trans_dict[trans_id][0] != sl[0]):
                print >> sys.stderr, 'WARNING: trans_id %s exists on chr %s\n ignore occurance on %s' % (trans_id, trans_dict[trans_id][0], sl[0])
            else:
                trans_dict[trans_id][2].append((int(sl[3]), int(sl[4])))
        except KeyError:
            trans_dict[trans_id] = [sl[0], sl[6],[(int(sl[3]), int(sl[4]))]]
    print >> log, "... parsed %i exons in %i transcripts\n" % (counter, len(trans_dict.keys())) 

    ### check annotation sanity
    ### sort exons in a transcript by start
    for key in trans_dict:
        trans_dict[key][2] = sorted(trans_dict[key][2], key = lambda x: x[0])

    ### parse fasta file to get genome sequence
    print >> log, "Parsing genome from %s ..." % options.fasta
    genome = dict()
    curr_chr = None
    for line in open(options.fasta,'r'):
        if line[0] == '>':
            if curr_chr:
                genome[curr_chr] = ''.join([x for x in genome[curr_chr]])
            curr_chr = options.chr_prefix + line[1:].strip()
            genome[curr_chr] = []
            print >> log, "... chr %s ... " % curr_chr
        else:
            genome[curr_chr].append(line.strip())
    if curr_chr:
        genome[curr_chr] = ''.join([x for x in genome[curr_chr]])
    print >> log, "... done\n"
            
    ### assemble transcriptome fasta
    try:
        fd_out = open(options.outfile, 'w')
        for trans_id in trans_dict:
            chr = trans_dict[trans_id][0]
            strand = trans_dict[trans_id][1]
            mrna = ''
            for exon in trans_dict[trans_id][2]:    
                mrna += genome[chr][exon[0] - 1:exon[1]]
            if strand == '-':   
                mrna = genome_utils.reverse_complement(mrna)

            print >> fd_out, ">%s" % trans_id 
            print >> fd_out, mrna
    except:
        pdb.set_trace()
    fd_out.close()


if __name__ == "__main__":
    main()
