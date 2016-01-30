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

#ifndef __UTILS_H__
#define __UTILS_H__

#include <string>
#include <cstdio>
#include <vector>
#include <deque>
#include <set>
#include <cmath>

#include "Alignment.h"
#include "GeneralData.h"
#include "config.h"

using namespace std;

extern GeneralData* genData;
extern Config* conf;

FILE* open_bam_pipe_in(std::string & in_fname);
FILE* open_bam_pipe_out(std::string & out_fname);

void parse_cigar(string cigar, vector<char> &operations, vector<int> &sizes);

double intron_penalty(vector<unsigned int> &intron_coverage);

double get_variance_global();

double get_variance(vector<vector<vector<unsigned long> > > &exon_coverage, vector<vector<set<unsigned long> > > &genome_pos);

vector<unsigned int> alter_coverage(vector<unsigned int> &source, unsigned int window_left, unsigned int window_right, bool is_positive);

string update_line_flag(char* line, bool is_best);

void parse_header(char* sl);

void get_single_loss(vector<Alignment>::iterator candidate, double &loss, bool debug = false);

void get_paired_loss(vector<Alignment>::iterator candidate_left, vector<Alignment>::iterator candidate_right, double &loss, bool debug = false);

bool compare_pair(vector<Alignment>::iterator candidate_left, vector<Alignment>::iterator candidate_right, vector<Alignment>::iterator best_left, vector<Alignment>::iterator best_right, double &loss, double &gain, bool debug = false);

bool compare_single(vector<Alignment>::iterator candidate, vector<Alignment>::iterator best, double &loss, double &gain, bool debug = false);

set<vector<Alignment>::iterator> filter_alignments(vector<Alignment> &aligns);

void  get_plifs_from_file();

double compute_mip_loss(double observed_cov, double predicted_cov, unsigned long segment_len = 0);

void prepare_mip_objective();

void add_zero_segments();

void compute_coverage_loss(vector<pair<vector<Alignment>::iterator,bool> > aligns, vector<vector<vector<unsigned long> > > &cov_keep, vector<vector<vector<unsigned long> > > &cov_change, vector<vector<set<unsigned long> > > &genome_pos);

bool pair_is_valid(vector<Alignment>::iterator first, vector<Alignment>::iterator second);

void parse_annotation();

void check_sorted_input(); 

#endif
