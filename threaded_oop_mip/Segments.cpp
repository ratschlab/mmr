
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
        char strand = '+';

        bool is_segment = true;

        while(sl != NULL) {
            // 0 -> s (segment) or i (intron)
            if (idx == 0) {
                is_segment = (strcmp("s", sl) == 0);
            } else if (idx == 1) { // 1 -> chr
               chr = genData->chr_num[string(sl)];
            } else if (idx == 2) { // 2 -> strand
                strand = *sl;
                if (strand == '.' && conf->strand_specific) {
                    fprintf(stderr, "Option for strand-specificity chosen, but given segments are unstranded!\n");
                    exit(2);
                } else if (strand == '-' && !conf->strand_specific) {
                    fprintf(stderr, "Given segment file seems strand specific, but option strand-specific is off. Please use '.' as strand in the segment file, if data is unstranded.\n");
                    exit(2);
                }
                if (strand == '.')
                    strand = '+';
            } else if (idx == 3) { // 3 -> start
                start = is_segment?atoi(sl):atoi(sl) + 1;
                //start-=1; // 1-based -> 0-based
            } else if (idx == 4) { // 4 -> stop
                stop = is_segment?atoi(sl):atoi(sl) - 1;
                //stop-=1; // 1-based -> 0-based
            } else if (idx == 5) { // 5 -> expected weight
                expectation = atof(sl);
            } else {
                break;
            }
            sl = strtok(NULL, "\t");
            idx++;
        }

        // construct entry in maps
        Segment* seg = new Segment(start, (unsigned int) (stop - start + 1), chr, strand, expectation);
        pair<unsigned char, unsigned char> chr_strand = pair<unsigned char, unsigned char>(chr, strand);
        if (is_segment) {
            if (this->exons.find(chr_strand) != this->exons.end()) {
                this->exons[chr_strand].insert(pair<long, long>(start, segment_id));
                this->exons[chr_strand].insert(pair<long, long>(stop, segment_id));
            } else {
                map<long, long> tmp_map;
                tmp_map.insert(pair<long, long>(start, segment_id));
                tmp_map.insert(pair<long, long>(stop, segment_id));
                this->exons.insert(pair<pair<unsigned char, unsigned char>, map<long, long> >(chr_strand, tmp_map));
            }
            this->exon_ids.insert(pair<long, Segment*>(segment_id, seg));
            segment_id++;
        } else {
            if (introns.find(chr_strand) != introns.end()) {
                this->introns[chr_strand].insert(pair<long, long>(start, intron_id));
                //this->introns[chr].insert(pair<long, long>(stop, intron_id));
            } else {
                multimap<unsigned long, unsigned long> tmp_map;
                tmp_map.insert(pair<unsigned long, unsigned long>(start, intron_id));
                //tmp_map.insert(pair<long, long>(stop, intron_id));
                this->introns.insert(pair<pair<unsigned char, unsigned char>, multimap<unsigned long, unsigned long> >(chr_strand, tmp_map));
            }
            this->introns_by_ids.insert(pair<long, Segment*>(intron_id, seg));
            intron_id++;
        }
    }
    if (conf->verbose)
        fprintf(stdout, "\t...parsed %i exon segments\n\t...parsed %i intron segments\n", (int) segment_id, (int) intron_id);
    fclose(infile);
}


pair<double, double> Segments::get_exon_segment_loss(vector<Alignment>::iterator alignment, set<unsigned long> overlap_region, bool debug) {

    // overlapping segments
    // pair of segment ID and length of overlap to the alignment (number of overlapping positions)
    map<unsigned long, set<unsigned int> > curr_exon_ids;
    vector<unsigned long> curr_intron_ids;

    // get all blocks of the current alignment
    // block coordinates are closed intervals!
    vector<pair<unsigned long, unsigned long> > exon_blocks;
    alignment->get_blocks(exon_blocks);

    // iterate over all available blocks and add the IDs of overlapping exon segments to the common list curr_exon_ids 
    vector<pair<unsigned long, unsigned long> >::iterator block;
    pair<unsigned char, unsigned char> chr_strand = pair<unsigned char, unsigned char>(alignment->chr, alignment->strand);
    for (block = exon_blocks.begin(); block != exon_blocks.end(); block++) {
        // lower is the first element in the ordered list of segment starts and stops
        // that is greater or equal to the block start
        map<long, long>::iterator lower = this->exons[chr_strand].lower_bound(block->first);
        // upper is the first element in the ordered list of segment starts and stops
        // that is greater to the block end
        map<long, long>::iterator upper = this->exons[chr_strand].upper_bound(block->second);
        if (upper != this->exons[chr_strand].end())
            upper++;
        
        // iterate over all overlapping segments
        multimap<long, long>::iterator i;
        set<long> handled_ids;
        for (i = lower; i != upper && i != this->exons[chr_strand].end(); i++) {
            //fprintf(stdout, "ali start: %i ali end: %i\n", alignment->start, alignment->get_end());
            //fprintf(stdout, "block start: %i block end: %i\n", block->first, block->second);
            //fprintf(stdout, "segment start: %i segment end: %i\n", this->exon_ids[i->second]->start, this->exon_ids[i->second]->start + this->exon_ids[i->second]->length - 1);
            //unsigned int affected_pos = min(alignment->get_end(), block->second) - max(alignment->start, block->first) + 1;
            
            //if (curr_exon_ids.find(i->second) == curr_exon_ids.end()) {
            if (handled_ids.find(i->second) == handled_ids.end()) {
                set<unsigned int> affected_pos;
                for (unsigned int v_i = max(this->exon_ids[i->second]->start, block->first); v_i <= min(this->exon_ids[i->second]->start + this->exon_ids[i->second]->length - 1, block->second); v_i++) {
                    affected_pos.insert(v_i);
                }
                curr_exon_ids.insert(pair<unsigned long, set<unsigned int>>(i->second, affected_pos));
                handled_ids.insert(i->second);
                //unsigned int affected_pos = min(this->exon_ids[i->second]->start + this->exon_ids[i->second]->length - 1, block->second) - max(this->exon_ids[i->second]->start, block->first) + 1;
                //curr_exon_ids.insert(pair<unsigned long, unsigned int>(i->second, affected_pos));
            }
        }
    }

    // get intron information from alignment blocks
    // intron IDs are inferred directly from the exon blocks
    if (exon_blocks.size() > 1) {
        //for (block = exon_blocks.begin(); block != exon_blocks.end(); block++) {
        //    fprintf(stdout, "block start: %i block end: %i\n", (unsigned int)block->first, (unsigned int) block->second);
        //}
        for (block = exon_blocks.begin() + 1; block != exon_blocks.end(); block++) {
            // get range of intron starts fitting to the previous block end
            pair<multimap<unsigned long, unsigned long>::iterator, multimap<unsigned long, unsigned long>::iterator> it_range = this->introns[chr_strand].equal_range((*(block-1)).second + 1);
            // check if any intron fits the gap between previous block and current block
            for (multimap<unsigned long, unsigned long>::iterator it = it_range.first; it != it_range.second; it++) {
                if (debug) {
                        fprintf(stdout, "checked intron from %lu to %lu\n", (*(block-1)).second + 1, block->first - 1);
                        fprintf(stdout, "with exp intron from %lu to %lu\n", this->introns_by_ids[it->second]->start, this->introns_by_ids[it->second]->start+this->introns_by_ids[it->second]->length - 1);
                    }
                    if ((*(block-1)).second + this->introns_by_ids[it->second]->length == block->first - 1) {
                        curr_intron_ids.push_back(it->second);
                        break;
                    }
                }
            }
        }

    // compute coverage and loss
    vector<unsigned int>::iterator cov_idx;

    // return value is >=  0 -> at least one overlapping segment exists or mip unpredicted regions are forced zero
    // return value is == -1 -> no overlapping segments exist
    // decision is made below
    double loss_with = 0.0;
    double loss_without = 0.0;

    if (curr_exon_ids.size() > 0) {
        // iterate over all exonic segments
        if (debug)
            fprintf(stdout, "Iterating over %i overlapping exonic segments:\n", (int) curr_exon_ids.size());
        map<unsigned long, set<unsigned int> >::iterator id;
        for (id = curr_exon_ids.begin(); id != curr_exon_ids.end(); id++) {
            unsigned long seg_start = this->exon_ids[id->first]->start;
            unsigned long seg_end = seg_start + this->exon_ids[id->first]->length;
            unsigned long seg_cov = 0;
            unsigned long overlap_len = 0;
            pair<unsigned char, unsigned char> cov_chr = pair<unsigned char, unsigned char>(alignment->chr, alignment->strand);
            for (cov_idx = genData->coverage_map[cov_chr].begin() + seg_start; cov_idx != genData->coverage_map[cov_chr].begin() + seg_end; cov_idx++) {
                // seg_cov stores the total coverage of the segment
                // overlap_cov stores the coverage of positions that come from the overlap between the possible multimapping locations
                //if (overlap_region.size() > 0 && overlap_region.find((cov_idx - genData->coverage_map[alignment->chr].begin()) + seg_start) != overlap_region.end())
                unsigned long ov_pos = cov_idx - genData->coverage_map[cov_chr].begin();
                if (overlap_region.size() > 0 && overlap_region.find(ov_pos) != overlap_region.end() && id->second.find(ov_pos) != id->second.end())
                    overlap_len += 1; 
                seg_cov += (unsigned long) (*cov_idx); 
            }
            if (debug)
                fprintf(stdout, "seg id: %lu seg start: %lu seg end: %lu seg strand: %c\n", id->first, seg_start, seg_end - 1, this->exon_ids[id->first]->strand);
            
           // compute mean segment coverage under consideration of the overlap
           // double mean_cov_with = 0.0;
           // double mean_cov_without = 0.0;
            unsigned long seg_cov_with = 0;
            unsigned long seg_cov_without = 0;
            //fprintf(stdout, "seg cov: %i affected pos: %i seg len: %i\n", seg_cov, id->second, this->exon_ids[id->first]->length);
            if (debug) 
                fprintf(stdout, "seg_cov: %lu affected pos: %i overlap_len: %lu\n", seg_cov, (int) id->second.size(), overlap_len);
            
            if (alignment->is_best) {
                seg_cov_with = seg_cov;
                seg_cov_without = seg_cov - id->second.size() + overlap_len;
            } else {
                seg_cov_with = seg_cov + id->second.size() - overlap_len;
                seg_cov_without = seg_cov;
            }
            if (debug)
                fprintf(stdout, "seg_cov_with: %lu seg_cov_without: %lu\n", seg_cov_with, seg_cov_without);

            // compute loss with and without current alignment
            loss_with += compute_mip_loss(seg_cov_with, this->exon_ids[id->first]->expectation, this->exon_ids[id->first]->length);
            loss_without += compute_mip_loss(seg_cov_without, this->exon_ids[id->first]->expectation, this->exon_ids[id->first]->length);
            //fprintf(stdout, "overlap cov: %i\n", overlap_cov);
            //fprintf(stdout, "exon_loss with %f; coverage with: %f; coverage expected: %f\n", loss_with, mean_cov_with, this->exon_ids[id->first]->expectation); 
            //fprintf(stdout, "exon_loss without %f; coverage without: %f; coverage expected: %f\n", loss_without, mean_cov_without, this->exon_ids[id->first]->expectation); 
        }
    }
    else if (!conf->use_mip_variance) {
        if (debug) {
            fprintf(stdout, "Did not find overlapping segments!\n");
            alignment->print();
        }
        loss_with = 0.0;
        loss_without = 0.0;
    } else {
        loss_with = -1.0;
        loss_without = -1.0;
    }

    if (curr_intron_ids.size() > 0) {
        // iterate over all intronic segments
        vector<unsigned long>::iterator id_int;
        if (debug)
            fprintf(stdout, "Iterating over %i overlapping exonic segments:\n", (int) curr_intron_ids.size());
        for (id_int = curr_intron_ids.begin(); id_int != curr_intron_ids.end(); id_int++) {
            unsigned int intron_cov = 0;
            unsigned long intron_start = this->introns_by_ids[*id_int]->start;
            unsigned long intron_end = intron_start + this->introns_by_ids[*id_int]->length - 1;
            map< pair<unsigned long, unsigned long>, unsigned int>::iterator it = genData->intron_coverage_map[pair<unsigned char, unsigned char>(alignment->chr, alignment->strand)].find(pair<unsigned long, unsigned long>(intron_start, intron_end));
            if (it != genData->intron_coverage_map[pair<unsigned char, unsigned char>(alignment->chr, alignment->strand)].end())
                intron_cov = it->second;
            if (debug)
                fprintf(stdout, "intron id: %lu intron start: %lu intron end: %lu intron strand: %c\n", *id_int, intron_start, intron_end - 1, this->introns_by_ids[*id_int]->strand);

            double intron_cov_with = 0.0;
            double intron_cov_without = 0.0;
            if (alignment->is_best) {
                intron_cov_with = intron_cov;
                intron_cov_without = intron_cov - 1;
            } else {
                intron_cov_with = intron_cov + 1;
                intron_cov_without = intron_cov;
            }
            // compute loss with and without current alignment
            loss_with += compute_mip_loss(intron_cov_with, this->introns_by_ids[*id_int]->expectation);
            loss_without += compute_mip_loss(intron_cov_without, this->introns_by_ids[*id_int]->expectation);
            //fprintf(stdout, "intron_loss with %f; coverage with: %f; coverage expected: %f\n", loss_with, intron_cov_with, this->introns_by_ids[*id_int]->expectation); 
            //fprintf(stdout, "intron_loss without %f; coverage without: %f; coverage expected: %f\n", loss_without, intron_cov_without, this->introns_by_ids[*id_int]->expectation); 
        }
    }

    pair<double, double> result (loss_with, loss_without);
    return result;
}

double Segments::get_total_loss() {
    
    double total_loss = 0.0;
    unsigned long total_ex_cov = 0;

   // map<long, double> loss_by_id;
  //  map<long, unsigned long> cov_by_id;

    // get total exon segment loss
    for (map<long, Segment*>::iterator it = this->exon_ids.begin(); it != this->exon_ids.end(); it++) {
        unsigned long curr_cov = 0;
        pair<unsigned char, unsigned char> cov_chr = pair<unsigned char, unsigned char>(it->second->chr, it->second->strand);
        for (vector<unsigned int>::iterator cov_idx = genData->coverage_map[cov_chr].begin() + it->second->start; cov_idx != genData->coverage_map[cov_chr].begin() + it->second->start + it->second->length; cov_idx++) {
            curr_cov += (unsigned long) (*cov_idx); 
        }
        //total_loss += compute_mip_loss(curr_cov / it->second->length, it->second->expectation);
        double curr_loss = compute_mip_loss(curr_cov, it->second->expectation, it->second->length);
        total_loss += curr_loss;
        total_ex_cov += curr_cov;

    //    loss_by_id.insert(pair<long, double>(it->first, curr_loss));
    //    cov_by_id.insert(pair<long, unsigned long>(it->first, curr_cov));
    }

    // get total intron segment loss
    for (map<long, Segment*>::iterator it = this->introns_by_ids.begin(); it != this->introns_by_ids.end(); it++) {
        pair<unsigned char, unsigned char> cov_chr = pair<unsigned char, unsigned char>(it->second->chr, it->second->strand);
        map< pair<unsigned long, unsigned long>, unsigned int>::iterator i_cov = genData->intron_coverage_map[cov_chr].find(pair<unsigned long, unsigned long>(it->second->start, it->second->start + it->second->length - 1));
        double curr_loss = 0.;
        if (i_cov != genData->intron_coverage_map[cov_chr].end())
            curr_loss = compute_mip_loss(i_cov->second, it->second->expectation);
        else
            curr_loss = compute_mip_loss(0, it->second->expectation);
        total_loss += curr_loss;
  //      loss_by_id.insert(pair<long, double>(it->first, curr_loss));
  //      cov_by_id.insert(pair<long, unsigned long>(it->first, i_cov != genData->intron_coverage_map[cov_chr].end()?i_cov->second:0));
    }

/*    // look for changes segments
    if (!genData->loss_by_segment.empty()) {
        map<long, double>::iterator it2 = genData->loss_by_segment.begin();
        map<long, unsigned long>::iterator it3 = cov_by_id.begin();
        for (map<long, double>::iterator it = loss_by_id.begin(); it != loss_by_id.end(); it++, it2++, it3++) {
            if (it->second != it2->second) {
                fprintf(stdout, "changed id: %li from: %f to %f\n", it->first, it2->second, it->second);
                fprintf(stdout, " \t start: %lu length: %i chr: %i strand:%c curr cov: %lu\n", this->exon_ids[it->first]->start, this->exon_ids[it->first]->length, this->exon_ids[it->first]->chr, this->exon_ids[it->first]->strand, it3->second); 
            }
        }
        fprintf(stdout, "\n");
    } 
    if (conf->num_threads < 2)
        genData->loss_by_segment = loss_by_id;

*/
    return total_loss;
}
