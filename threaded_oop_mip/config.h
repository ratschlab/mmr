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

        bool verbose;
        bool print_best_only;
        bool use_pair_info;
        bool use_variants;
        bool pre_filter;
        bool parse_complete;
        bool use_mip_objective;
        
        unsigned int window_size;
        unsigned int iterations;
        unsigned int num_threads;
        unsigned int max_fifo_size;
        unsigned char filter_distance;

        int intron_offset;

        string infile;
        string outfile;
        string segments;

        double insert_size;
        double insert_dev;
    private:
        void init();
};
#endif
