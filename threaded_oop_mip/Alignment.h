#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include <cstring>
#include <string>
#include <set>
#include "GeneralData.h"

using namespace std;

class Alignment {

public:
    unsigned char chr;
    unsigned long start;
    vector<char> operations;
    vector<int> sizes;
    bool is_best;
    unsigned char edit_ops;
    unsigned char quality;
    bool reversed;

    Alignment(const unsigned char a, const unsigned long b, const vector<char> c, const vector<int> cc, const bool d, const unsigned char e, const unsigned char f, bool g) :
        chr(a), start(b), operations(c), sizes(cc), is_best(d), edit_ops(e) , quality(f), reversed(g) {}

    Alignment() :
        chr(0), start(0), is_best(false), edit_ops(0), quality(0), reversed(false) {}

    ~Alignment() {};

    bool operator==(const Alignment &other);

    string fill(char* sl, unsigned char &pair);

    void get_coverage(unsigned int window_size, vector<unsigned short> &intron_cov, vector<unsigned short> &exon_cov);
    void update_coverage_map(bool positive);
    unsigned int get_end();
    
    bool comparator(const Alignment &left, const Alignment &right); 

    bool is_spliced();

    vector< pair<unsigned long, unsigned int> > get_intron_coords();

    void get_blocks(vector<pair<unsigned long, unsigned long> > &blocks);

    set<unsigned long> get_overlap(Alignment &other);
};
#endif
