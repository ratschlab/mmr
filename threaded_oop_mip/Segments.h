#ifndef __SEGMENTS_H__
#define __SEGMENTS_H__

using namespace std;

#include <map>
#include <set>
#include <utility>
#include <string>

#include "Alignment.h"
#include "GeneralData.h"
#include "Segment.h"

class Segments {

public:

    // stores exon start or stap as key, segment id as value
    multimap<long, long> exons;
    // stores segment id as key, pair of segment start and len as value
    map<long, Segment> exon_ids;
    // stores intron start or stap as key, segment id as value
    multimap<long, long> introns;
    // stores segment id as key, pair of segment start and len as value
    map<long, Segment> intron_ids;

    Segments() {};
    ~Segments() {};

    void get_from_file(string &filename);

    pair<float, float> get_exon_segment_loss(vector<Alignment>::iterator alignment, set<unsigned long> overlap_region, bool is_best);
    pair<float, float> get_intron_segment_loss(vector<Alignment>::iterator alignment);
};
#endif
