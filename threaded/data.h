#ifndef __DATA_H__
#define __DATA_H__

#include "alignment.h"

#include <string>

using namespace std;

struct GeneralData {

    // general coverage structures
    map <string, int> chr_num;
    vector <unsigned int> chr_size;
    map <int, vector<unsigned short> > coverage_map;

    // best hit maps (will stay empty for batch setting)
    unordered_map <string, size_t, hash<string> > best_left;
    unordered_map <string, size_t, hash<string> > best_right;

};

struct ThreadData_Batch {

    // reads
    unordered_map <string, vector<alignment>, hash<string> > read_map_left;
    unordered_map <string, vector<alignment>, hash<string> > read_map_right;

    // active sets
    unordered_map <string, vector<vector<alignment>::iterator>, hash<string> > active_left_pair;
    unordered_map <string, vector<vector<alignment>::iterator>, hash<string> > active_right_pair;
    list<unordered_map <string, vector<alignment> >::iterator > active_left_single;
    list<unordered_map <string, vector<alignment> >::iterator > active_right_single;

};

struct ThreadData_Online {

    vector<alignment> left_reads;
    vector<alignment> right_reads;
    
    string last_id;
};

#endif
