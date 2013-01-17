#ifndef __SEGMENT_H__
#define __SEGMENT_H__

using namespace std;

class Segment {

    public:

    unsigned long start;
    unsigned int length;
    unsigned char chr;
    unsigned char strand;
    float expectation;

    Segment(const unsigned long a, const unsigned int b, const unsigned char c,  const unsigned char d, const float e) :
        start(a), length(b), chr(c), strand(d), expectation(e) {}

    Segment():
        start(0), length(0), chr(1), strand('+'), expectation(0.0) {}

    ~Segment() {};
    
    bool operator==(const Segment &other) {
        return (start == other.start && length == other.length && chr == other.chr && strand == other.strand);
    };

};
#endif
