
#include <algorithm>

#include "Segments.h"
#include "Utils.h"

extern GeneralData* genData;

void Segments::get_from_file(string & filename) {
    FILE* infile = fopen(filename.c_str(), "r");
    char* ret;
    char line[1000];

    while (true) {
        ret = fgets(line, sizeof(line), infile);

        if (!ret)
            break;

        if (line[0] == '#')
            continue;
    }

    char* sl = strtok(line, "\t");

    // TODO More code here
    int idx = 0;
    long segment_id = 0;
    long intron_id = 0;
    unsigned int chr = 1;
    long start = 0;
    long stop = 0;
    float expectation = 0.0;

    bool is_segment = true;

    while(sl != NULL) {
        string tmp_sl = string(sl); 
        // 0 -> s (segment) or i (intron)
        if (idx == 0) {
            is_segment = (strcmp("s", sl) == 0);
        } else if (idx == 1) { // 1 -> chr
           chr = genData->chr_num[tmp_sl];
        } else if (idx == 2) { // 2 -> strand
            // ignored for now  TODO
        } else if (idx == 3) { // 3 -> start
            start = atoi(sl);
        } else if (idx == 4) { // 4 -> stop
            stop = atoi(sl);
        } else if (idx == 5) { // 5 -> expected weight
            expectation = atof(sl);
        } else {
            break;
        }
        idx++;
    }

    // construct entry in maps
    Segment* seg = new Segment(start, (unsigned int) (stop - start + 1), chr, expectation);
    if (is_segment) {
        if (this->exons.find(chr) != this->exons.end()) {
            this->exons[chr].insert(pair<long, long>(start, segment_id));
            this->exons[chr].insert(pair<long, long>(stop, segment_id));
        } else {
            multimap<long, long> tmp_map;
            tmp_map.insert(pair<long, long>(start, segment_id));
            tmp_map.insert(pair<long, long>(stop, segment_id));
            this->exons.insert(pair<unsigned int, multimap<long, long> >(chr, tmp_map));
        }
        this->exon_ids.insert(pair<long, Segment*>(segment_id, seg));
        segment_id++;
    } else {
        if (introns.find(chr) != introns.end()) {
            this->introns[chr].insert(pair<long, long>(start, intron_id));
            this->introns[chr].insert(pair<long, long>(stop, intron_id));
        } else {
            multimap<long, long> tmp_map;
            tmp_map.insert(pair<long, long>(start, intron_id));
            tmp_map.insert(pair<long, long>(stop, intron_id));
            this->introns.insert(pair<unsigned int, multimap<long, long> >(chr, tmp_map));
        }
        this->intron_ids.insert(pair<long, Segment*>(intron_id, seg));
        intron_id++;
    }
}


pair<double, double> Segments::get_exon_segment_loss(vector<Alignment>::iterator alignment, set<unsigned long> overlap_region, bool is_best) {

    // overlapping segments
    // pair of segment ID and length of overlap to the alignment (number of overlapping positions)
    map<unsigned long, unsigned int> curr_ids;

    // get all blocks of the current alignment
    vector<pair<unsigned long, unsigned long> > blocks;
    alignment->get_blocks(blocks);

    // iterate over all available blocks and add the IDs of overlapping exon segments to the common list curr_ids 
    vector<pair<unsigned long, unsigned long> >::iterator block;
    for (block = blocks.begin(); block != blocks.end(); block++) {
        // lower is the first element in the ordered list of segment starts and stops
        // that is greater or equal to the block start
        multimap<long, long>::iterator lower = this->exons[alignment->chr].lower_bound(block->first);
        multimap<long, long>::iterator upper = this->exons[alignment->chr].lower_bound(block->second);
        // upper is the last element in the range of values that correspond to the first key 
        // in the list of ordered segment starts and stops that is greater or equal to the alignment end 
        upper = this->exons[alignment->chr].equal_range(upper->first).second;

        // iterate over all overlapping segments
        multimap<long, long>::iterator i;
        for (i = lower; i != upper; i++) {
            unsigned int affected_pos = min(alignment->get_end(), block->second) - max(alignment->start, block->first) + 1;
            if (curr_ids.find(i->second) != curr_ids.end())
                curr_ids[i->second] += affected_pos;
            else
                curr_ids.insert(pair<unsigned long, unsigned int>(i->second, affected_pos));
        }
    }

    // compute coverage and loss
    vector<unsigned short>::iterator cov_idx;

    // return value is >=  0 -> at least one overlapping segment exists
    // return value is == -1 -> no overlapping segments exist
    double loss_with = 0.0;
    double loss_without = 0.0;
    if (curr_ids.size() <= 0) {
        loss_with = -1.0;
        loss_without = -1.0;
    }

    // iterate over all segments
    map<unsigned long, unsigned int>::iterator id;
    for (id = curr_ids.begin(); id != curr_ids.end(); id++) {
        unsigned long seg_start = this->exon_ids[id->first]->start;
        unsigned long seg_end = seg_start + this->exon_ids[id->first]->length;
        unsigned int seg_cov = 0;
        unsigned int overlap_cov = 0;
        for (cov_idx = genData->coverage_map[alignment->chr].begin() + seg_start; cov_idx != genData->coverage_map[alignment->chr].begin() + seg_end; cov_idx++) {
            // seg_cov stores the total coverage of the segment
            // overlap_cov stores the coverage of positions that come from the overlap between the possible multimapping locations
            if (overlap_region.size() > 0 && overlap_region.find((cov_idx - genData->coverage_map[alignment->chr].begin()) + seg_start) != overlap_region.end())
                overlap_cov += (unsigned int) (*cov_idx);
            seg_cov += (unsigned int) (*cov_idx); 
        }
        
        // compute mean segment coverage under consideration of the overlap
        double mean_cov_with = 0.0;
        double mean_cov_without = 0.0;
        if (is_best) {
            mean_cov_with = (double) seg_cov / this->exon_ids[id->first]->length;
            mean_cov_without = (double) (seg_cov - id->second + overlap_cov) / this->exon_ids[id->first]->length;
        } else {
            mean_cov_with = (double) (seg_cov + id->second - overlap_cov) / this->exon_ids[id->first]->length;
            mean_cov_without = (double)seg_cov / this->exon_ids[id->first]->length;
        }

        // compute loss with and without current alignment
        loss_with += compute_loss(mean_cov_with, this->exon_ids[id->first]->expectation);
        loss_without += compute_loss(mean_cov_without, this->exon_ids[id->first]->expectation);
    }

    pair<float, float> result (loss_with, loss_without);
    return result;
}

set<long> Segments::get_affected_intron_segs(vector<Alignment>::iterator alignment) {
    
    set<long> result;

    // check, if alignment has at least one intron
    if (! alignment->is_spliced())
        return result;

    // get list of introns
    vector< pair< unsigned long, unsigned int> > align_introns = alignment->get_intron_coords();

    // iterate over all introns
    vector< pair< unsigned long, unsigned int> >::iterator it;
    for (it = align_introns.begin(); it != align_introns.end(); it++) {
        multimap<long, long>::iterator start = this->introns[alignment->chr].find(it->first);
        multimap<long, long>::iterator stop = this->introns[alignment->chr].find(it->first + (unsigned long) it->second - 1);
        
        if (start != this->introns[alignment->chr].end() && stop != this->introns[alignment->chr].end() && start == stop) {
            result.insert(start->second);
        }
    }

    return result;
}
