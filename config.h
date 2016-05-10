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
        bool print_unmapped;
        bool use_pair_info;
        bool use_variants;
        bool pre_filter;
        bool parse_complete;
        bool use_mip_objective;
        bool use_mip_variance;
        bool strand_specific;
        bool zero_unpred;
        bool burn_in;
        bool debug;
        bool fast_mutex;
        bool take_non_secondary_only; 
        bool use_brkpts;
        bool no_sort_check;
        
        unsigned int window_size;
        unsigned int iterations;
        unsigned int num_threads;
        unsigned int max_fifo_size;
        unsigned int read_len;
        unsigned int max_list_length;
        unsigned int max_pair_list_length;
        unsigned char filter_distance;
        unsigned char trim_id;

        int intron_offset;
        int max_gen_frag_size;

        string infile;
        string outfile;
        string segmentfile;
        string annotation;
        string lossfile;
        string samtools;

        double insert_dev;

        int iteration;
        double last_loss;

    private:
        void init();
};
#endif
