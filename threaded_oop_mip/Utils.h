#ifndef __UTILS_H__
#define __UTILS_H__

#include <string>
#include <cstdio>
#include <vector>
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

double get_variance(vector<unsigned int> &exon_coverage, vector<unsigned int> &intron_coverage);

vector<unsigned int> alter_coverage(vector<unsigned int> &source, unsigned int window_left, unsigned int window_right, bool is_positive);

string update_line_flag(char* line, bool is_best);

void parse_header(char* sl);

bool compare_pair(vector<Alignment>::iterator candidate_left, vector<Alignment>::iterator candidate_right, vector<Alignment>::iterator best_left, vector<Alignment>::iterator best_right, double &loss, bool debug = false);

bool compare_single(vector<Alignment>::iterator candidate, vector<Alignment>::iterator best, double &loss, bool debug = false);

set<vector<Alignment>::iterator> filter_alignments(vector<Alignment> &aligns);

void  get_plifs_from_file();

double compute_mip_loss(double observed_cov, double predicted_cov, unsigned long segment_len = 0);

void prepare_mip_objective();

void add_zero_segments();

#endif
