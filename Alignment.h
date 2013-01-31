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

    void fill_coverage_vectors(vector<unsigned long> &cov_keep, vector<unsigned long> &cov_change, set<unsigned int> &genome_pos, unsigned long first_start, bool is_curr_best);

    void fill_coverage_vector(vector<unsigned long> &cov_keep);

    void alter_coverage_vector(vector<vector<unsigned long> > &cov_change, vector<set<unsigned long> > &genome_pos, bool is_curr_best);

    void update_coverage_map(bool positive);

    unsigned long get_end();
    
    bool compare_edit_ops(const Alignment &left, const Alignment &right); 

    bool is_spliced();

    vector< pair<unsigned long, unsigned int> > get_intron_coords();

    void get_blocks(vector<pair<unsigned long, unsigned long> > &blocks);

    set<unsigned long> get_overlap(Alignment &other);

    set<unsigned long> get_genome_pos(unsigned int window_size = 0);

    unsigned int get_length();

    void determine_gaps();

    void clear();
    
    void print();
};
#endif
