#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <cstring>
#include <string>

using namespace std;

class Config {
    public:
        Config();
        Config(int argc, char *argv[]);

        ~Config();
       
        void print_usage(string prog_name);
        void print_call(string prog_name);

        bool verbose;
        bool print_best_only;
        bool use_pair_info;
        bool use_variants;
        bool pre_filter;
        bool parse_complete;
        bool use_mip_objective;
        bool use_mip_variance;
        bool strand_specific;
        
        unsigned int window_size;
        unsigned int iterations;
        unsigned int num_threads;
        unsigned int max_fifo_size;
        unsigned int read_len;
        unsigned char filter_distance;

        int intron_offset;

        string infile;
        string outfile;
        string segmentfile;
        string lossfile;

        double insert_size;
        double insert_dev;

        int iteration;
        double last_loss;

    private:
        void init();
};
#endif
