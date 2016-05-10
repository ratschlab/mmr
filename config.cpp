/*

  This file is part of MMR, the Read Multi-Mapper Resolution tool.

  MMR is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This software is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  A copy of the GNU General Public License is distributed with 
  MMR (file LICENSE).  If not, see <http://www.gnu.org/licenses/>.

  Written 2010-2015 by 
    Andre Kahles <akahles@cbio.mskcc.org>
    Jonas Behr <jonas_behr@web.de>
    Gunnar R\"atsch <raetsch@cbio.mskcc.org>
 
  This work was funded by the Max Planck Society,
  the German Research Foundation (DFG RA1894/2-1)
  and Memorial Sloan Kettering Cancer Center.

*/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "config.h"

// default constructor
Config::Config() {
    init();
}

Config::Config(int argc, char *argv[]) {

    init();
    int i = 0;
    if (argc == 1) {
        print_usage(string(argv[0]));
        exit(-1);
    }
    while (i < argc) {
        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
            verbose = true;
        } else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--best-only")) {
            print_best_only = true;
        } else if (!strcmp(argv[i], "-u") || !strcmp(argv[i], "--keep-unmapped")) {
            print_unmapped = true;
        } else if (!strcmp(argv[i], "-V") || !strcmp(argv[i], "--use-variants")) {
            use_variants = true;
        } else if (!strcmp(argv[i], "-C") || !strcmp(argv[i], "--init-secondary")) {
            take_non_secondary_only = false;
        } else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--pair-usage")) {
            use_pair_info = true;
        } else if (!strcmp(argv[i], "-S") || !strcmp(argv[i], "--strand-specific")) {
            strand_specific = true;
        } else if (!strcmp(argv[i], "-m") || !strcmp(argv[i], "--mitie-objective")) {
            use_mip_objective = true;
        } else if (!strcmp(argv[i], "-M") || !strcmp(argv[i], "--mitie-variance")) {
            use_mip_variance = true;
        } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--max-fragment-size")) {
            max_gen_frag_size = atoi(argv[++i]);
//        } else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--insert-dev")) {
//            insert_dev = (double) atof(argv[++i]);
        } else if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "--pre-filter-off")) {
            pre_filter = false;
        } else if (!strcmp(argv[i], "-P") || !strcmp(argv[i], "--parse-complete")) {
            parse_complete = true;
        } else if (!strcmp(argv[i], "-F") || !strcmp(argv[i], "--filter-dist")) {
            filter_distance = (unsigned char) atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--trim-id")) {
            trim_id = (unsigned char) atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-R") || !strcmp(argv[i], "--read-len")) {
            read_len = (unsigned char) atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-w") || !strcmp(argv[i], "--windowsize")) {
            window_size = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-a") || !strcmp(argv[i], "--annotation")) {
            annotation = std::string(argv[++i]);
            use_brkpts = true;
        } else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--threads")) {
            num_threads = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-L") || !strcmp(argv[i], "--max-list-length")) {
            max_list_length = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-A") || !strcmp(argv[i], "--max-pair-list-length")) {
            max_pair_list_length = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-I") || !strcmp(argv[i], "--iterations")) {
            iterations = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--debug")) {
            debug = true;
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(std::string(argv[0]));
            exit(0);
        } else if (!strcmp(argv[i], "-o")) {
            outfile = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-T") || !strcmp(argv[i], "--path-to-samtools")) {
            samtools = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-z") || !strcmp(argv[i], "--zero-expect-unpred")) {
            zero_unpred = true;
        } else if (!strcmp(argv[i], "-B") || !strcmp(argv[i], "--burn-in")) {
            burn_in = true;
        } else if (!strcmp(argv[i], "--no-sort-check")) {
            no_sort_check = true;
        } else if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "--segmentfile")) {
            segmentfile = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--lossfile")) {
            lossfile = std::string(argv[++i]);
        } else {
            if (argv[i][0] == '-') {
                fprintf(stderr, "\nERROR: Unknown option %s\n", argv[i]);
                print_usage(std::string(argv[0]));
                exit(-1);
            }
            else {
                infile = string(argv[i]);
            }
        }
        i++;
    }
}

// nothing to destroy
Config::~Config() {}

void Config::print_usage(std::string prog_name) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "\t%s -o OUTFILE [options] IN_BAM\n\n", prog_name.c_str());
    fprintf(stderr, "Available Options:\n\n");
    // Input and threads
    fprintf(stderr, "\n\tInput handling and paralellization:\n");
    fprintf(stderr, "\t-P --parse-complete \tparse complete file into memory [off]\n");
    fprintf(stderr, "\t-t --threads \t\tnumber of threads to use (must be > 2) [1]\n");
    fprintf(stderr, "\t-S --strand-specific \talignments are strand specific [off]\n");
    fprintf(stderr, "\t-C --init-secondary  \tchoose initial alignment also from secondary lines (flag 256) [off]\n");
    fprintf(stderr, "\t   --no-sort-check  \tinput files are not checked to be sorted by read ID [off]\n");
    // Filter options
    fprintf(stderr, "\n\tInput file filtering:\n");
    fprintf(stderr, "\t-f --pre-filter-off \tswitch off pre filter for alignments that have F more edit ops than the best [on]\n");
    fprintf(stderr, "\t-F --filter-dist [INT]\tfilter distance F for pre-filter [1]\n");
    fprintf(stderr, "\t-V --use-variants \tuse variant alignments for filtering (different edit op count,\n");
    fprintf(stderr, "\t\t\t\trequires XG and XM Tag in alignment files) [off]\n");
    fprintf(stderr, "\t-L --max-list-length [INT]\tmax length of alignment list per read (after filtering) [1000]\n");
    fprintf(stderr, "\t-r --trim-id \t\ttrim this many positions from the end of each read ID [0]\n");
    // Paired Alignment options
    fprintf(stderr, "\n\tPaired alignment handling:\n");
    fprintf(stderr, "\t-p --pair-usage \tpre use pair information in the reads [off]\n");
//    fprintf(stderr, "\t-d --insert-dev \tallowed deviation from insert size (times insert size) [0.4]\n");
    fprintf(stderr, "\t-i --max-fragment-size \tupper limit of GENOMIC fragment length [1 000 000]\n");
    fprintf(stderr, "\t-A --max-pair-list-length [INT]\tmax no of valid pairs before not using pair modus [10000]\n");
    // Output Options
    fprintf(stderr, "\n\tOutput handling:\n");
    fprintf(stderr, "\t-b --best-only \t\tprint only best alignment [off]\n");
    fprintf(stderr, "\t-u --keep-unmapped \tprint unmapped reads from input [off]\n");
    // options for variance optimization
    fprintf(stderr, "\n\tOptions for using the variance optimization:\n");
    fprintf(stderr, "\t-w --windowsize  [INT]\tsize of coverage window around read [20]\n");
    fprintf(stderr, "\t-I --iterations  [INT]\tnumber of iterations to smooth the coverage [5]\n");
    fprintf(stderr, "\t-a --annotation \tannotation file in GTF format to infer segment boundaries []\n");
    fprintf(stderr, "\t-B --burn-in \t\tuse the first iteration to fill coverage map only [off]\n");
    // MIP Options
    fprintf(stderr, "\n\tOptions for using the MiTie objective for smoothing:\n");
    fprintf(stderr, "\t-m --mitie-objective \tuse objective from MiTie instead of local variance [off]\n");
    fprintf(stderr, "\t-s --segmentfile \tMiTie segment file required for MiTie optimization []\n");
    fprintf(stderr, "\t-l --lossfile \t\tMiTie loss parameter file required for MiTie optimization []\n");
    fprintf(stderr, "\t-R --read-len  [INT]\taverage length of the reads [75]\n");
    fprintf(stderr, "\t-M --mitie-variance \tuse variance smoothing for regions with no MiTie prediction [off]\n");
    fprintf(stderr, "\t-z --zero-expect-unpred \tinitializes all covered but not predicted positions with expectation 0.0 [off]\n");
    // General options
    fprintf(stderr, "\n\tGeneral:\n");
    fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
    fprintf(stderr, "\t-h --help \t\tprint usage info\n");
    fprintf(stderr, "\n");
}

void Config::print_call(std::string prog_name) {
    fprintf(stdout, "%s has been started with the following parameters:\n\n", prog_name.c_str());
    fprintf(stdout, "\t input file:           %s\n", infile.c_str()); 
    fprintf(stdout, "\t output file:          %s\n", outfile.c_str()); 
    fprintf(stdout, "\t strand specific:      %s\n", strand_specific?"yes":"no");
    fprintf(stdout, "\t pre filter:           %s\n", pre_filter?"on":"off");
    fprintf(stdout, "\t init on secondary:    %s\n", take_non_secondary_only?"off":"on");
    if (pre_filter) {
        fprintf(stdout, "\t filter dist:          %i\n", filter_distance);
        fprintf(stdout, "\t use variants:         %s\n", use_variants?"on":"off");
    }
    fprintf(stdout, "\t max list length:      %i\n", max_list_length);
    fprintf(stdout, "\t pair usage:           %s\n", use_pair_info?"on":"off");
    if (use_pair_info) {
        fprintf(stdout, "\t max frag size size:   %i\n", max_gen_frag_size);
        fprintf(stdout, "\t max pair list length: %i\n", max_pair_list_length);
//        fprintf(stdout, "\t insert size std dev:  %.2f\n", insert_dev);
    }
    if (trim_id > 0)
        fprintf(stdout, "\t trim read id by:      %i\n", trim_id);
    fprintf(stdout, "\t print best only:      %s\n", print_best_only?"on":"off");
    fprintf(stdout, "\t print unmapped:       %s\n", print_unmapped?"on":"off");
    fprintf(stdout, "\t iterations:           %i\n", iterations);
    fprintf(stdout, "\t 1 iteration burn in:  %s\n", burn_in?"on":"off");
    if (use_brkpts) 
        fprintf(stdout, "\t annotation file:      %s\n", annotation.c_str()); 
    fprintf(stdout, "\t threads:              %i\n", num_threads);
    if (use_mip_variance || ! use_mip_objective) {
        fprintf(stdout, "\t window size:          %i\n", window_size);
    } else {
        fprintf(stdout, "\t use MiTie objective:    %s\n", use_mip_objective?"on":"off");
        fprintf(stdout, "\t segment file:         %s\n", segmentfile.c_str()); 
        fprintf(stdout, "\t loss file:            %s\n", lossfile.c_str()); 
        fprintf(stdout, "\t read length:          %i\n", read_len); 
        fprintf(stdout, "\t zero for unpred seg:  %s\n", zero_unpred?"yes":"no");
        fprintf(stdout, "\t use variance if no MiTie-segment overlaps: %s\n", use_mip_variance?"on":"off");
    }
}

// PRIVATE
void Config::init() {
    verbose = false;
    print_best_only = false;
    print_unmapped = false;
    use_pair_info = false;
    pre_filter = true;
    parse_complete = false;
    use_mip_objective = false;
    strand_specific = false;
    debug = false;
    fast_mutex = true;
    window_size = 20;
    iterations = 5;
    filter_distance = 1;
    trim_id = 0;
    intron_offset = 5;
    max_gen_frag_size = 1000000;
//    insert_dev = 0.4;
    num_threads = 1;
    max_fifo_size = 5000;
    max_list_length = 1000;
    max_pair_list_length = 10000;
    use_variants = false;
    use_mip_variance = false;
    segmentfile = string();
    annotation = string();
    lossfile = string();
    read_len = 75;
    zero_unpred = false;
    use_brkpts = false;
    burn_in = false;
    no_sort_check = false;
    samtools = "samtools";
    take_non_secondary_only = true;

    iteration = 0;
    last_loss = 0.0;
}

