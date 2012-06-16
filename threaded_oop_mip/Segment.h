#ifndef __SEGMENT_H__
#define __SEGMENT_H__

using namespace std;

class Segment {

    public:

    unsigned long start;
    unsigned int length;
    float expectation;

    Segment(const unsigned long a, const unsigned int b, const float c) :
        start(a), length(b), expectation(c) {}

    Segment():
        start(0), length(0), expectation(0.0) {}

    ~Segment() {};
    
    bool operator==(const Segment &other) {
        return (start == other.start && length == other.length);
    };

};
#endif
