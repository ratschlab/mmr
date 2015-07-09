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

#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include <string>
#include <vector>
#include <set>
#include <map>

using namespace std;

class Alignment {

public:
    unsigned int chr;
    unsigned long start;
    vector<char> operations;
    vector<int> sizes;
    bool is_best;
    unsigned char edit_ops;
    unsigned char quality;
    unsigned char strand;
    bool reversed;
    bool is_secondary;
    set<unsigned int> deletions;
    map<unsigned int, int> insertions;

    Alignment(const unsigned int a, const unsigned long b, const vector<char> c, const vector<int> cc, const bool d, const unsigned char e, const unsigned char f, const unsigned char g, const bool h, const bool i, const set<unsigned char> deletions, const map<unsigned int, int> insertions) :
        chr(a), start(b), operations(c), sizes(cc), is_best(d), edit_ops(e) , quality(f), strand(g), reversed(h), is_secondary(i) {}

    Alignment() :
        chr(0), start(0), is_best(false), edit_ops(0), quality(0), strand('+'), reversed(false), is_secondary(false) {}

    ~Alignment() {};

    bool operator==(const Alignment &other);

    string fill(char* sl, unsigned char &pair, bool &unmapped);

    void fill_coverage_vector(vector<vector<unsigned long> > &cov_keep);

    void alter_coverage_vector(vector<vector<vector<unsigned long> > > &cov_change, vector<vector<set<unsigned long> > > &genome_pos, bool is_curr_best);

    void update_coverage_map(bool positive);

    unsigned long get_end();
    
    bool compare_edit_ops(const Alignment &left, const Alignment &right); 

    bool is_spliced();

    vector< pair<unsigned long, unsigned int> > get_intron_coords();

    void get_blocks(vector<pair<unsigned long, unsigned long> > &blocks);

    set<unsigned long> get_overlap(Alignment &other);

    vector<set<unsigned long> > get_genome_pos(unsigned int window_size = 0);

    set<unsigned long> get_exon_pos(unsigned int window_size = 0);

    unsigned int get_length();

    void determine_gaps();

    void clear();
    
    void print();
};
#endif
