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
