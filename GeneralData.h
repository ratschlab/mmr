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

#ifndef __GENERALDATA_H__
#define __GENERALDATA_H__

#include <string>
#include <map>
#include <unordered_map>
#include <vector>

#include "Segments.h"

using namespace std;

struct GeneralData {

    // general coverage structures
    map <string, unsigned int> chr_num;
    vector <unsigned int> chr_size;
    vector <unsigned int> chr_size_cum;
    map <pair<unsigned int, unsigned char>, vector<unsigned int> > coverage_map; // contains coverage vector for each chr/strand pair
    map <pair<unsigned int, unsigned char>, vector<bool> > breakpoint_map; // contains breakpoint coordinates induced by annotated segment boundaries

    // intron coverage map (will only be used for Mitie optimization)
    // outer map: keys -> chr/strand pair  values -> intron map
    // inner map: keys -> start/end pair (closed interval)  values -> coverages
    map <pair<unsigned int, unsigned char>, map< pair<unsigned long, unsigned long>, unsigned int> > intron_coverage_map;

    // best hit maps (will stay empty for batch setting)
    unordered_map <string, size_t, hash<string> > best_left;
    unordered_map <string, size_t, hash<string> > best_right;

    // plifs to compute loss in case of mip objective
    map< double, vector<double> > plifs;

    // data structure containing segment graph information
    Segments segments;

    // count data 
    double total_loss;
    double total_gain;
    unsigned int num_altered;
    unsigned int num_tested;

    // loss by segment id
    map<long, double> loss_by_segment;

};
#endif

