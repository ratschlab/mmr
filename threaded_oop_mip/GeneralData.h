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
    map <string, int> chr_num;
    vector <unsigned int> chr_size;
    map <int, vector<unsigned short> > coverage_map;

    // intron coverage map (will only be used for MIP optimization)
    map <int, map< pair<unsigned long, unsigned long>, unsigned int> > intron_coverage_map;

    // best hit maps (will stay empty for batch setting)
    unordered_map <string, size_t, hash<string> > best_left;
    unordered_map <string, size_t, hash<string> > best_right;

    // plifs to compute loss in case of mip objective
    map< double, vector<double> > plifs;

    // data structure containing segment graph information
    Segments segments;

    // count data 
    double total_loss;
    unsigned int num_altered;

};
#endif

