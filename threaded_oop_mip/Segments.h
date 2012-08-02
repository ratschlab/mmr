#ifndef __SEGMENTS_H__
#define __SEGMENTS_H__

using namespace std;

#include <map>
#include <set>
#include <utility>
#include <string>

#include "Alignment.h"
#include "Segment.h"

class Segments {

public:


    // outer map: chr as key, inner map as value
    // inner map: stores exon start or stop as key, segment id as value
    map<unsigned int, map<long, long> > exons;
    // stores segment id as key, segment object as value
    map<long, Segment*> exon_ids;
    // outer map: chr as key, inner map as value
    // inner map: stores intron start as key, segment id as value
    map<unsigned int, multimap<unsigned long, unsigned long> > introns;
    // stores segment id as key, segment object as value
    map<long, Segment*> introns_by_ids;

    Segments() {};
    ~Segments() {};

    void get_from_file();

    pair<double, double> get_exon_segment_loss(vector<Alignment>::iterator alignment, set<unsigned long> overlap_region, bool is_best);
    pair<double, double> get_intron_segment_loss(vector<Alignment>::iterator alignment);

    double get_total_loss();
};
#endif
