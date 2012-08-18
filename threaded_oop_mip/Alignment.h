#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include <cstring>
#include <string>
#include <vector>
#include <set>

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
    unsigned char strand;
    bool reversed;

    Alignment(const unsigned char a, const unsigned long b, const vector<char> c, const vector<int> cc, const bool d, const unsigned char e, const unsigned char f, const unsigned char g, bool h) :
        chr(a), start(b), operations(c), sizes(cc), is_best(d), edit_ops(e) , quality(f), strand(g), reversed(h) {}

    Alignment() :
        chr(0), start(0), is_best(false), edit_ops(0), quality(0), strand('+'), reversed(false) {}

    ~Alignment() {};

    bool operator==(const Alignment &other);

    string fill(char* sl, unsigned char &pair);

    //pair<double, double> get_variance_loss(set<unsigned long> overlap_region);
    pair<double, double> get_variance_loss(set<unsigned long> covered_pos, set<unsigned long> not_covered_pos);

    void update_coverage_map(bool positive);
    unsigned long get_end();
    
    bool comparator(const Alignment &left, const Alignment &right); 

    bool is_spliced();

    vector< pair<unsigned long, unsigned int> > get_intron_coords();

    void get_blocks(vector<pair<unsigned long, unsigned long> > &blocks);

    set<unsigned long> get_overlap(Alignment &other);

    set<unsigned long> get_genome_pos();

    void clear();
    
    void print();
};
#endif
