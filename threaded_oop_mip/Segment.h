#ifndef __SEGMENT_H__
#define __SEGMENT_H__

using namespace std;

class Segment {

    public:

    unsigned long start;
    unsigned int length;
    unsigned int chr;
    float expectation;

    Segment(const unsigned long a, const unsigned int b, const unsigned int c, const float d) :
        start(a), length(b), chr(c), expectation(d) {}

    Segment():
        start(0), length(0), chr(1), expectation(0.0) {}

    ~Segment() {};
    
    bool operator==(const Segment &other) {
        return (start == other.start && length == other.length && chr == other.chr);
    };

};
#endif
