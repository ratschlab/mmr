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

#ifndef __BATCHDATA_H__
#define __BATCHDATA_H__

#include <vector>
#include <map>
#include <unordered_map>
#include <list>

#include "Utils.h"
#include "GeneralData.h"
#include "Alignment.h"

using namespace std;

class BatchData {

public:
    // reads
    unordered_map <string, vector<Alignment>, hash<string> > read_map_left;
    unordered_map <string, vector<Alignment>, hash<string> > read_map_right;

    // active sets
    unordered_map <string, vector<vector<Alignment>::iterator>, hash<string> > active_left_pair;
    unordered_map <string, vector<vector<Alignment>::iterator>, hash<string> > active_right_pair;
    list<unordered_map <string, vector<Alignment> >::iterator > active_left_single;
    list<unordered_map <string, vector<Alignment> >::iterator > active_right_single;

    // constructor 
    BatchData() {};

    // destructor 
    ~BatchData() {};

    void parse_file();
    void pre_filter_alignment_maps();
    double get_total_min_loss();
    void get_active_read_set();

    unsigned int smooth_coverage_map_single(unsigned int &num_ambiguous);
    unsigned int smooth_coverage_map_paired(unsigned int &num_ambiguous);

private:
    unsigned int smooth_coverage_map_single_wrapper(list<unordered_map <string, vector<Alignment> >::iterator > &active_reads, unsigned int &num_ambiguous);
};
#endif
