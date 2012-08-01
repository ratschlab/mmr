
#include <algorithm>

#include "Segments.h"
#include "Utils.h"

extern GeneralData* genData;
extern Config* conf;

void Segments::get_from_file() {
    FILE* infile = fopen(conf->segmentfile.c_str(), "r");
    char* ret;
    char line[1000];
    long segment_id = 0;
    long intron_id = 0;

    while (true) {
        ret = fgets(line, sizeof(line), infile);

        if (!ret)
            break;

        if (line[0] == '#')
            continue;

        char* sl = strtok(line, "\t");

        int idx = 0;
        unsigned int chr = 1;
        long start = 0;
        long stop = 0;
        float expectation = 0.0;

        bool is_segment = true;

        while(sl != NULL) {
            // 0 -> s (segment) or i (intron)
            if (idx == 0) {
                is_segment = (strcmp("s", sl) == 0);
            } else if (idx == 1) { // 1 -> chr
               chr = genData->chr_num[string(sl)];
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
            sl = strtok(NULL, "\t");
            idx++;
        }

        // construct entry in maps
        Segment* seg = new Segment(start, (unsigned int) (stop - start + 1), chr, expectation);
        if (is_segment) {
            if (this->exons.find(chr) != this->exons.end()) {
                this->exons[chr].insert(pair<long, long>(start, segment_id));
                this->exons[chr].insert(pair<long, long>(stop, segment_id));
            } else {
                map<long, long> tmp_map;
                tmp_map.insert(pair<long, long>(start, segment_id));
                tmp_map.insert(pair<long, long>(stop, segment_id));
                this->exons.insert(pair<unsigned int, map<long, long> >(chr, tmp_map));
            }
            this->exon_ids.insert(pair<long, Segment*>(segment_id, seg));
            segment_id++;
        } else {
            if (introns.find(chr) != introns.end()) {
                this->introns[chr].insert(pair<long, long>(start, intron_id));
                //this->introns[chr].insert(pair<long, long>(stop, intron_id));
            } else {
                multimap<unsigned long, unsigned long> tmp_map;
                tmp_map.insert(pair<unsigned long, unsigned long>(start, intron_id));
                //tmp_map.insert(pair<long, long>(stop, intron_id));
                this->introns.insert(pair<unsigned int, multimap<unsigned long, unsigned long> >(chr, tmp_map));
            }
            this->introns_by_ids.insert(pair<long, Segment*>(intron_id, seg));
            intron_id++;
        }
    }
    if (conf->verbose)
        fprintf(stdout, "\t...parsed %i exon segments\n\t...parsed %i intron segments\n", (int) segment_id, (int) intron_id);
    fclose(infile);
}


pair<double, double> Segments::get_exon_segment_loss(vector<Alignment>::iterator alignment, set<unsigned long> overlap_region, bool is_best) {

    // overlapping segments
    // pair of segment ID and length of overlap to the alignment (number of overlapping positions)
    map<unsigned long, unsigned int> curr_exon_ids;
    vector<unsigned long> curr_intron_ids;

    // get all blocks of the current alignment
    // block coordinates are closed intervals!
    vector<pair<unsigned long, unsigned long> > exon_blocks;
    alignment->get_blocks(exon_blocks);

    // iterate over all available blocks and add the IDs of overlapping exon segments to the common list curr_exon_ids 
    vector<pair<unsigned long, unsigned long> >::iterator block;
    for (block = exon_blocks.begin(); block != exon_blocks.end(); block++) {
        // lower is the first element in the ordered list of segment starts and stops
        // that is greater or equal to the block start
        map<long, long>::iterator lower = this->exons[alignment->chr].lower_bound(block->first);
        // upper is the first element in the ordered list of segment starts and stops
        // that is greater to the block end
        map<long, long>::iterator upper = this->exons[alignment->chr].upper_bound(block->second);
        
        // iterate over all overlapping segments
        multimap<long, long>::iterator i;
        for (i = lower; i != upper; i++) {
            unsigned int affected_pos = min(alignment->get_end(), block->second) - max(alignment->start, block->first) + 1;
            if (curr_exon_ids.find(i->second) != curr_exon_ids.end())
                curr_exon_ids[i->second] += affected_pos;
            else
                curr_exon_ids.insert(pair<unsigned long, unsigned int>(i->second, affected_pos));
        }
    }

    // get intron information from alignment blocks
    // intron IDs are inferred directly from the exon blocks
    if (exon_blocks.size() > 1) {
        for (block = exon_blocks.begin() + 1; block != exon_blocks.end(); block++) {
            // get range of intron starts fitting to the previous block end
            pair<multimap<unsigned long, unsigned long>::iterator, multimap<unsigned long, unsigned long>::iterator> it_range = this->introns[alignment->chr].equal_range((*(block-1)).second + 1);
            for (multimap<unsigned long, unsigned long>::iterator it = it_range.first; it != it_range.second; it++) {
                // check if any intron fits the gap between previous block and current block
                if ((*(block-1)).second + this->introns_by_ids[it->second]->length == block->first - 1) {
                    curr_intron_ids.push_back(it->second);
                    break;
                }
            }
        }
    }

    // compute coverage and loss
    vector<unsigned short>::iterator cov_idx;

    // return value is >=  0 -> at least one overlapping segment exists or mip unpredicted regions are forced zero
    // return value is == -1 -> no overlapping segments exist
    // decision is made below
    double loss_with = 0.0;
    double loss_without = 0.0;

    if (curr_exon_ids.size() > 0) {
        // iterate over all exonic segments
        map<unsigned long, unsigned int>::iterator id;
        for (id = curr_exon_ids.begin(); id != curr_exon_ids.end(); id++) {
            unsigned long seg_start = this->exon_ids[id->first]->start;
            unsigned long seg_end = seg_start + this->exon_ids[id->first]->length;
            unsigned int seg_cov = 0;
            unsigned int overlap_cov = 0;
            for (cov_idx = genData->coverage_map[alignment->chr].begin() + seg_start; cov_idx != genData->coverage_map[alignment->chr].begin() + seg_end; cov_idx++) {
                // seg_cov stores the total coverage of the segment
                // overlap_cov stores the coverage of positions that come from the overlap between the possible multimapping locations
                if (overlap_region.size() > 0 && overlap_region.find((cov_idx - genData->coverage_map[alignment->chr].begin()) + seg_start) != overlap_region.end())
                    overlap_cov += 1; //(unsigned int) (*cov_idx);
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
            loss_with += compute_mip_loss(mean_cov_with, this->exon_ids[id->first]->expectation);
            loss_without += compute_mip_loss(mean_cov_without, this->exon_ids[id->first]->expectation);
        }
    }
    else if (!conf.use_mip_variance) {
        // iterate over blocks instead of segments
        for (block = exon_blocks.begin(); block != exon_blocks.end(); block++) {
            unsigned int block_cov = 0;
            unsigned int overlap_cov = 0;
            for (cov_idx = genData->coverage_map[alignment->chr].begin() + block->first; cov_idx != genData->coverage_map[alignment->chr].begin() + block->second + 1; cov_idx++) {
                // block_cov stores the total coverage of all blocks
                // overlap_cov stores the coverage of positions that come from the overlap between the possible multimapping locations
                if (overlap_region.size() > 0 && overlap_region.find((cov_idx - genData->coverage_map[alignment->chr].begin()) + block->first) != overlap_region.end())
                    overlap_cov += 1;
                block_cov += (unsigned int) (*cov_idx); 
            }
            // compute mean block coverage under consideration of the overlap
            double mean_cov_with = 0.0;
            double mean_cov_without = 0.0;
            unsigned long block_len = block->second - block->first + 1;
            if (is_best) {
                mean_cov_with = (double) block_cov / block_len;
                mean_cov_without = (double) (block_cov - block_len + overlap_cov) / block_len;
            } else {
                mean_cov_with = (double) (block_cov + block_len - overlap_cov) / block_len;
                mean_cov_without = (double) block_cov / block_len;
            }

            // compute loss with and without current alignment
            loss_with += compute_mip_loss(mean_cov_with, 0);
            loss_without += compute_mip_loss(mean_cov_without, 0);
    } else {
        loss_with = -1.0;
        loss_without = -1.0;
    }

    if (curr_intron_ids.size() > 0) {
        // iterate over all intronic segments
        vector<unsigned long>::iterator id_int;
        for (id_int = curr_intron_ids.begin(); id_int != curr_intron_ids.end(); id_int++) {
            unsigned int intron_cov = 0;
            unsigned long intron_start = this->introns_by_ids[*id_int]->start;
            unsigned long intron_end = intron_start + this->introns_by_ids[*id_int]->length - 1;
            map< pair<unsigned long, unsigned long>, unsigned int>::iterator it = genData->intron_coverage_map[alignment->chr].find(pair<unsigned long, unsigned long>(intron_start, intron_end));
            if (it != genData->intron_coverage_map[alignment->chr].end())
                intron_cov = it->second;

            double intron_cov_with = 0.0;
            double intron_cov_without = 0.0;
            if (is_best) {
                intron_cov_with = intron_cov;
                intron_cov_without = intron_cov - 1;
            } else {
                intron_cov_with = intron_cov + 1;
                intron_cov_without = intron_cov;
            }
            // compute loss with and without current alignment
            loss_with += compute_mip_loss(intron_cov_with, this->introns_by_ids[*id_int]->expectation);
            loss_without += compute_mip_loss(intron_cov_without, this->introns_by_ids[*id_int]->expectation);
        }
    }
    else if (!conf.use_mip_variance) {
        // iterate over 
        if (exon_blocks.size() > 1) {
            for (block = exon_blocks.begin() + 1; block != exon_blocks.end(); block++) {
                unsigned int intron_cov = 0;
                map< pair<unsigned long, unsigned long>, unsigned int>::iterator it = genData->intron_coverage_map[alignment->chr].find(pair<unsigned long, unsigned long>((*(block-1)).second + 1, block->first - 1));
                if (it != genData->intron_coverage_map[alignment->chr].end())
                    intron_cov = it->second;
                
                double intron_cov_with = 0.0;
                double intron_cov_without = 0.0;
                if (is_best) {
                    intron_cov_with = intron_cov;
                    intron_cov_without = intron_cov - 1;
                } else {
                    intron_cov_with = intron_cov + 1;
                    intron_cov_without = intron_cov;
                }
                // compute loss with and without current alignment
                loss_with += compute_mip_loss(intron_cov_with, 0);
                loss_without += compute_mip_loss(intron_cov_without, 0);
            }
        }
    }

    pair<float, float> result (loss_with, loss_without);
    return result;
}
