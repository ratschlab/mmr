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
