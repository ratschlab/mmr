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
        } else if (!strcmp(argv[i], "-b")) {
            print_best_only = true;
        } else if (!strcmp(argv[i], "-V") || !strcmp(argv[i], "--use_variants")) {
            use_variants = true;
        } else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--pair-usage")) {
            use_pair_info = true;
        } else if (!strcmp(argv[i], "-m") || !strcmp(argv[i], "--mip-objective")) {
            use_mip_objective = true;
        } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--insert-size")) {
            insert_size = (double) atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--insert-dev")) {
            insert_size = (double) atof(argv[++i]);
        } else if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "--pre-filter")) {
            pre_filter = true;
        } else if (!strcmp(argv[i], "-P") || !strcmp(argv[i], "--parse-complete")) {
            parse_complete = true;
        } else if (!strcmp(argv[i], "-F") || !strcmp(argv[i], "--filter-dist")) {
            filter_distance = (unsigned char) atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-w") || !strcmp(argv[i], "--windowsize")) {
            window_size = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--threads")) {
            num_threads = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-I") || !strcmp(argv[i], "--iterations")) {
            iterations = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(std::string(argv[0]));
            exit(0);
        } else if (!strcmp(argv[i], "-o")) {
            outfile = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "--segments")) {
            segments = std::string(argv[++i]);
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
    fprintf(stderr, "Available Options:\n");
    fprintf(stderr, "\t-w --windowsize  [INT]\tsize of coverage window around read [20]\n");
    fprintf(stderr, "\t-I --iterations  [INT]\tnumber of iterations to smooth the coverage [5]\n");
    fprintf(stderr, "\t-F --filter-dist [INT]\tfilter distance F for pre-filter [3]\n");
    fprintf(stderr, "\t-f --pre-filter \tpre filter all alignments that have F more edit ops than the best [off]\n");
    fprintf(stderr, "\t-m --mip-objective \tuse objective from MIP instead of local variance [off]\n");
    fprintf(stderr, "\t-s --segments \tsegment file required for mip optimization []\n");
    fprintf(stderr, "\t-P --parse-complete \tparse complete file into memory [off]\n");
    fprintf(stderr, "\t-p --pair-usage \tpre use pair information in the reads [off]\n");
    fprintf(stderr, "\t-i --insert-size \testimted insert size for paired end reads [200]\n");
    fprintf(stderr, "\t-t --threads \tnumber of threads to use (must be > 2) [1]\n");
    fprintf(stderr, "\t-d --insert-dev \tallowed deviation from insert size (times insert size) [0.4]\n");
    fprintf(stderr, "\t-b --best-only \t\tprint only best alignment [off]\n");
    fprintf(stderr, "\t-V --use-variants \t\tuse variant alignments (different edit op count,\n");
    fprintf(stderr, "\t\t\t\t\trequires XG and XM Tag in alignment files) [off]\n");
    fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
    fprintf(stderr, "\t-h --help \t\tprint usage info\n");
}


// PRIVATE
void Config::init() {
    verbose = false;
    print_best_only = false;
    use_pair_info = false;
    pre_filter = false;
    parse_complete = false;
    use_mip_objective = false;
    window_size = 20;
    iterations = 5;
    filter_distance = 3;
    intron_offset = 5;
    insert_size = 200.0;
    insert_dev = 0.4;
    num_threads = 1;
    max_fifo_size = 1000;
    use_variants = false;
    segments = string();
}

