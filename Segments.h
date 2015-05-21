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

#ifndef __SEGMENTS_H__
#define __SEGMENTS_H__

using namespace std;

#include <map>
#include <set>
#include <utility>
#include <string>

#include "Alignment.h"
#include "Segment.h"

class Segments {

public:


    // outer map: chr/strand pair as key, inner map as value
    // inner map: stores exon start or stop as key, segment id as value
    map<pair<unsigned char, unsigned char>, map<long, long> > exons;
    // stores segment id as key, segment object as value
    map<long, Segment*> exon_ids;
    // outer map: chr/strand pair as key, inner map as value
    // inner map: stores intron start as key, segment id as value
    map<pair<unsigned char, unsigned char>, multimap<unsigned long, unsigned long> > introns;
    // stores segment id as key, segment object as value
    map<unsigned long, Segment*> introns_by_ids;

    Segments() {};
    ~Segments() {
        for (map<long, Segment*>::iterator it = exon_ids.begin(); it != exon_ids.end(); it++)
            delete it->second;
        for (map<unsigned long, Segment*>::iterator it = introns_by_ids.begin(); it != introns_by_ids.end(); it++)
            delete it->second;
    };

    void get_from_file();

    pair<double, double> get_exon_segment_loss(vector<vector<Alignment>::iterator> alignments, set<unsigned long> invariant_pos, bool debug = false);

    double get_total_loss();
};
#endif
