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
    map <string, unsigned char> chr_num;
    vector <unsigned int> chr_size;
    map <pair<unsigned char, unsigned char>, vector<unsigned int> > coverage_map; // contains coverage vector for each chr/strand pair

    // intron coverage map (will only be used for Mitie optimization)
    // outer map: keys -> chr/strand pair  values -> intron map
    // inner map: keys -> start/end pair (closed interval)  values -> coverages
    map <pair<unsigned char, unsigned char>, map< pair<unsigned long, unsigned long>, unsigned int> > intron_coverage_map;

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

    // loss by segment id
    map<long, double> loss_by_segment;

};
#endif

