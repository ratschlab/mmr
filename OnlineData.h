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

#ifndef __ONLINEDATA_H__
#define __ONLINEDATA_H__

#include <string>
#include <vector>
#include <set>
#include <cstdio>
#include <time.h>

#include "GeneralData.h"
#include "Alignment.h"

using namespace std;

class OnlineData {

public:

    vector<Alignment> left_reads;
    vector<Alignment> right_reads;
    
    string last_id;

    OnlineData() {};

    ~OnlineData() {};

    void process_data_online(GeneralData* genData);

    void get_active_reads(string read_id, set<vector<Alignment>::iterator> &ignore_reads_left, set<vector<Alignment>::iterator> &ignore_reads_right, vector<vector<Alignment>::iterator> &active_left_reads, vector<vector<Alignment>::iterator> &active_right_reads, GeneralData* genData, bool &found_pairs);

    char* parse_file(FILE* infile, char* last_line, GeneralData* genData, unsigned int &counter, clock_t &start_clock, clock_t &start_time);
};
#endif
