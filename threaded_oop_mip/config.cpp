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
        } else if (!strcmp(argv[i], "-S") || !strcmp(argv[i], "--strand-specific")) {
            strand_specific = true;
        } else if (!strcmp(argv[i], "-m") || !strcmp(argv[i], "--mip-objective")) {
            use_mip_objective = true;
        } else if (!strcmp(argv[i], "-M") || !strcmp(argv[i], "--mip-variance")) {
            use_mip_variance = true;
        } else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--read-len")) {
            read_len = (double) atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--insert-size")) {
            insert_size = (double) atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--insert-dev")) {
            insert_dev = (double) atof(argv[++i]);
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
        } else if (!strcmp(argv[i], "-z") || !strcmp(argv[i], "--zero-expect-unpred")) {
            zero_unpred = true;
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
    // Filter options
    fprintf(stderr, "\n\tInput file filtering:\n");
    fprintf(stderr, "\t-f --pre-filter \tpre filter all alignments that have F more edit ops than the best [on]\n");
    fprintf(stderr, "\t-F --filter-dist [INT]\tfilter distance F for pre-filter [1]\n");
    fprintf(stderr, "\t-V --use-variants \tuse variant alignments for filtering (different edit op count,\n");
    fprintf(stderr, "\t\t\t\trequires XG and XM Tag in alignment files) [off]\n");
    // Paired Alignment options
    fprintf(stderr, "\n\tPaired alignment handling:\n");
    fprintf(stderr, "\t-p --pair-usage \tpre use pair information in the reads [off]\n");
    fprintf(stderr, "\t-d --insert-dev \tallowed deviation from insert size (times insert size) [0.4]\n");
    fprintf(stderr, "\t-i --insert-size \testimted insert size for paired end reads [200]\n");
    // Output Options
    fprintf(stderr, "\n\tOutput handling:\n");
    fprintf(stderr, "\t-b --best-only \t\tprint only best alignment [off]\n");
    // options for variance optimization
    fprintf(stderr, "\n\tOptions for using the variance optimization:\n");
    fprintf(stderr, "\t-w --windowsize  [INT]\tsize of coverage window around read [20]\n");
    fprintf(stderr, "\t-I --iterations  [INT]\tnumber of iterations to smooth the coverage [5]\n");
    // MIP Options
    fprintf(stderr, "\n\tOptions for using the MiTie objective for smoothing:\n");
    fprintf(stderr, "\t-m --mip-objective \tuse objective from MiTie instead of local variance [off]\n");
    fprintf(stderr, "\t-s --segmentfile \tsegment file required for mip optimization []\n");
    fprintf(stderr, "\t-l --lossfile \t\tloss parameter file required for mip optimization []\n");
    fprintf(stderr, "\t-r --read-len  [INT]\taverage length of the reads [75]\n");
    fprintf(stderr, "\t-M --mip-variance \tuse variance smoothing for regions with no MiTie prediction [off]\n");
    fprintf(stderr, "\t-z --zero-expect-unpred \tinitializes all covered but not predicted positions with expectation 0.0 [off]\n");
    // General options
    fprintf(stderr, "\n\tGeneral:\n");
    fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
    fprintf(stderr, "\t-h --help \t\tprint usage info\n");
}

void Config::print_call(std::string prog_name) {
    fprintf(stdout, "%s has been started with the following parameters:\n\n", prog_name.c_str());
    fprintf(stdout, "\t input file:           %s\n", infile.c_str()); 
    fprintf(stdout, "\t output file:          %s\n", outfile.c_str()); 
    fprintf(stdout, "\t strand specific:      %s\n", pre_filter?"yes":"no");
    fprintf(stdout, "\t pre filter:           %s\n", pre_filter?"on":"off");
    if (pre_filter) {
        fprintf(stdout, "\t filter dist:          %i\n", filter_distance);
        fprintf(stdout, "\t use variants:         %s\n", use_variants?"on":"off");
    }
    fprintf(stdout, "\t pair usage:           %s\n", use_pair_info?"on":"off");
    if (use_pair_info) {
        fprintf(stdout, "\t insert size:          %.2f\n", insert_size);
        fprintf(stdout, "\t insert size std dev:  %.2f\n", insert_dev);
    }
    fprintf(stdout, "\t print best only:      %s\n", print_best_only?"on":"off");
    if (use_mip_variance || ! use_mip_objective) {
        fprintf(stdout, "\t window size:          %i\n", window_size);
        fprintf(stdout, "\t iterations:           %i\n", iterations);
    } else {
        fprintf(stdout, "\t use mip objective:    %s\n", use_mip_objective?"on":"off");
        fprintf(stdout, "\t segment file:         %s\n", segmentfile.c_str()); 
        fprintf(stdout, "\t loss file:            %s\n", lossfile.c_str()); 
        fprintf(stdout, "\t read length:          %i\n", read_len); 
        fprintf(stdout, "\t zero for unpred seg:  %s\n", zero_unpred?"yes":"no");
        fprintf(stdout, "\t use variance if no mip-segment overlaps: %s\n", use_mip_variance?"on":"off");
    }
}

// PRIVATE
void Config::init() {
    verbose = false;
    print_best_only = false;
    use_pair_info = false;
    pre_filter = true;
    parse_complete = false;
    use_mip_objective = false;
    strand_specific = false;
    window_size = 20;
    iterations = 5;
    filter_distance = 1;
    intron_offset = 5;
    insert_size = 200.0;
    insert_dev = 0.4;
    num_threads = 1;
    max_fifo_size = 1000;
    use_variants = false;
    use_mip_variance = false;
    segmentfile = string();
    lossfile = string();
    read_len = 75;
    zero_unpred = false;

    iteration = 0;
    last_loss = 0.0;
}

