#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include <cstring>
#include <string>

using namespace std;

struct alignment {
    unsigned char chr;
    unsigned long start;
    string cigar;
    bool is_best;
    unsigned char edit_ops;
    unsigned char quality;
    bool reversed;

    alignment(const unsigned char a, const unsigned long b, const string c, const bool d, const unsigned char e, const unsigned char f, bool g) :
        chr(a), start(b), cigar(c), is_best(d), edit_ops(e) , quality(f), reversed(g) {}

    alignment() :
        chr(0), start(0), cigar(string("")), is_best(false), edit_ops(0), quality(0), reversed(false) {}

    bool operator==(const alignment &other) {
        return (chr == other.chr && start == other.start && cigar == other.cigar && reversed == other.reversed);
    }
};

bool alignment_comparator(const alignment &left, const alignment &right) {
        return left.edit_ops < right.edit_ops;
}
#endif
