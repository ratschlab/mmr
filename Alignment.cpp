/*

  This file is part of MMR, the Read Multi-Mapper Resolution tool.

  MMR is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This software is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  A copy of the GNU General Public License is distributed with 
  MMR (file LICENSE).  If not, see <http://www.gnu.org/licenses/>.

  Written 2010-2015 by 
    Andre Kahles <akahles@cbio.mskcc.org>
    Jonas Behr <jonas_behr@web.de>
    Gunnar R\"atsch <raetsch@cbio.mskcc.org>
 
  This work was funded by the Max Planck Society,
  the German Research Foundation (DFG RA1894/2-1)
  and Memorial Sloan Kettering Cancer Center.

*/

#include <algorithm>

#include "assert.h"
#include "Alignment.h"
#include "Utils.h"

extern GeneralData* genData;
extern Config* conf;

extern pthread_mutex_t mutex_coverage;
extern pthread_mutex_t mutex_int_coverage;
extern pthread_mutex_t mutex_fifo;
extern pthread_mutex_t mutex_done;
extern pthread_mutex_t mutex_best_left;
extern pthread_mutex_t mutex_best_right;

bool Alignment::operator==(const Alignment &other) {
    return (chr == other.chr && start == other.start && operations == other.operations && sizes == other.sizes && reversed == other.reversed);
}

string Alignment::fill(char* sl, unsigned char &pair, bool &unmapped) {

    int idx = 0;
    string id;

    while (sl != NULL) {
        if (idx == 0) { 
            id = sl;
            if (conf->trim_id > 0)
                id = id.substr(0, id.size() - conf->trim_id);
        } else if (idx == 1) {
            pair = (atoi(sl) & 128);
            this->reversed = ((atoi(sl) & 16) == 16);
            this->is_secondary = ((atoi(sl) & 256) == 256);
        } else if (idx == 2) {
            unmapped = (sl[0] == '*');
            if (unmapped)
                break;
            if (genData->chr_num.find(sl) == genData->chr_num.end()) {
                fprintf(stderr, "ERROR: Contig name not in header!\n Contig: %s\n\n", sl);
                exit(-1);
            } else {
                this->chr = genData->chr_num[sl];
            }
        } else if (idx == 3) {
            this->start = strtoul(sl, NULL, 0);
            this->start--; // 1-based --> 0-based
        } else if (idx == 4) {
            this->quality = (unsigned  char) atoi(sl);
        } else if (idx == 5) {
            parse_cigar(sl, this->operations, this->sizes);
        } else if (!conf->use_variants && strlen(sl) > 3 && sl[0] == 'N' && sl[1] == 'M' && sl[2] == ':') {
            this->edit_ops = atoi(sl+5); 
            if (!conf->strand_specific)
                break;
        } else if (conf->use_variants && strlen(sl) > 3 && sl[0] == 'X' && (sl[1] == 'M' || sl[1] == 'G') && sl[2] == ':') {
           this-> edit_ops += atoi(sl+5);
            if (!conf->strand_specific)
                break;
        } else if (conf->strand_specific && strlen(sl) > 3 && sl[0] == 'X' && sl[1] == 'S' && sl[2] == ':') {
            this->strand = *(sl+5);
        }
        sl = strtok(NULL, "\t");
        idx ++;
    }
    if (idx < 6) {
        return(string(""));
    }
    this->determine_gaps();
    return id;
}

void Alignment::fill_coverage_vector(vector<vector<unsigned long> > &cov_keep) {

    pair<unsigned int, unsigned char> chr_strand = make_pair(this->chr, this->strand);
    vector<unsigned int>::iterator cov_idx = genData->coverage_map[ chr_strand ].begin();
    vector<bool>::iterator brk_idx = genData->breakpoint_map[ chr_strand ].begin();

    size_t chrm_pos = 0;
    if (conf->window_size < this->start) {
        advance(cov_idx, this->start - conf->window_size);
        advance(brk_idx, this->start - conf->window_size);
        chrm_pos = (this->start - conf->window_size);
    }
    unsigned int offset = (conf->window_size >= this->start)?this->start:conf->window_size;

    // lock coverage map
    if (conf->fast_mutex)
        pthread_mutex_lock(&mutex_coverage);

    if (conf->debug) {
        this->print();
        fprintf(stderr, "cov keep before fill:\n");
        for (vector<vector<unsigned long> >::iterator u = cov_keep.begin(); u != cov_keep.end(); u++) {
            for (vector<unsigned long>::iterator v = u->begin(); v != u->end(); v++) {
                fprintf(stderr, "%2lu ", *v);
            }
        }
        fprintf(stderr, "\n");
    }


    // get coverage from preceding windows
    vector<unsigned long> curr_cov_segment;
    for (size_t i = 0; i < offset; i++) {
        if (cov_idx < genData->coverage_map[ chr_strand ].end()) {
            curr_cov_segment.push_back(*cov_idx);
            cov_idx++;
            chrm_pos++;
            brk_idx++;
            if (conf->use_brkpts && brk_idx != genData->breakpoint_map[ chr_strand ].end() && *brk_idx) {
                curr_cov_segment.clear();
            }
        } else {
            break;
        }
    }
    if (conf->fast_mutex)
        pthread_mutex_unlock(&mutex_coverage);

    // get coverage from alignment exons
    for (size_t i = 0; i < this->sizes.size(); i++) {
        switch (this->operations.at(i)) {
            case 'M': case 'D': case 'N': {   //curr_pos = distance(genome_pos.begin(), genome_pos.find(chrm_pos + genData->chr_size_cum.at(this->chr - 1)));
                                    if (conf->fast_mutex)
                                        pthread_mutex_lock(&mutex_coverage);
                                    if (this->operations.at(i) == 'N' && this->sizes.at(i) > (2 * (int) conf->window_size)) {
                                        size_t pos2brkpt = 0;
                                        for (size_t j = 0; j < conf->window_size; j++) {
                                            if (cov_idx < genData->coverage_map[ chr_strand ].end()) {
                                                curr_cov_segment.push_back(*cov_idx);
                                                cov_idx++;
                                                chrm_pos++;
                                                brk_idx++;
                                                pos2brkpt++;
                                                if (conf->use_brkpts && brk_idx != genData->breakpoint_map[ chr_strand ].end() && *brk_idx) {
                                                    cov_keep.push_back(curr_cov_segment);
                                                    curr_cov_segment.clear();
                                                    break;
                                                }
                                            }
                                        }
                                        cov_idx += (this->sizes.at(i) - conf->window_size - pos2brkpt);
                                        brk_idx += (this->sizes.at(i) - conf->window_size - pos2brkpt);
                                        chrm_pos += (this->sizes.at(i) - conf->window_size - pos2brkpt);
                                        if (curr_cov_segment.size() > 0) {
                                            cov_keep.push_back(curr_cov_segment);
                                            curr_cov_segment.clear();
                                        }
                                        for (size_t j = 0; j < conf->window_size; j++) {
                                            if (cov_idx < genData->coverage_map[ chr_strand ].end()) {
                                                curr_cov_segment.push_back(*cov_idx);
                                                cov_idx++;
                                                brk_idx++;
                                                chrm_pos++;
                                                if (conf->use_brkpts && brk_idx != genData->breakpoint_map[ chr_strand ].end() && *brk_idx) {
                                                    curr_cov_segment.clear();
                                                }
                                            }
                                        }
                                    } else {
                                        for (int j = 0; j < this->sizes.at(i); j++) {
                                            if (cov_idx < genData->coverage_map[ chr_strand ].end()) {
                                                curr_cov_segment.push_back(*cov_idx);
                                                cov_idx++;
                                                brk_idx++;
                                                chrm_pos++;
                                                if (conf->use_brkpts && brk_idx != genData->breakpoint_map[ chr_strand ].end() && *brk_idx) {
                                                    cov_keep.push_back(curr_cov_segment);
                                                    curr_cov_segment.clear();
                                                }
                                            }
                                        }
                                    };
                                    if (conf->fast_mutex)
                                        pthread_mutex_unlock(&mutex_coverage);
                                    break;
                                }
        }
    }

    // get coverage from following windows
    if (conf->fast_mutex)
        pthread_mutex_lock(&mutex_coverage);
    for (size_t i = 0; i < conf->window_size; i++) {
        if (cov_idx < genData->coverage_map[ chr_strand ].end()) {
            curr_cov_segment.push_back(*cov_idx);
            cov_idx++;
            brk_idx++;
            //if (conf->use_brkpts && cov_idx < genData->coverage_map[ chr_strand ].end() && genData->breakpoint_map[ chr_strand ].at(distance(cov_beg, cov_idx))) {
            if (conf->use_brkpts && brk_idx != genData->breakpoint_map[ chr_strand ].end() && *brk_idx) {
                cov_keep.push_back(curr_cov_segment);
                curr_cov_segment.clear();
                break;
            }
        }
        else {
            break;
        }
    }

    if (curr_cov_segment.size() > 0)
        cov_keep.push_back(curr_cov_segment);

    if (conf->debug) {
        fprintf(stderr, "cov keep after fill:\n");
        for (vector<vector<unsigned long> >::iterator v = cov_keep.begin(); v != cov_keep.end(); v++) {
            for (vector<unsigned long>::iterator u = v->begin(); u != v->end(); u++) {
                fprintf(stderr, "%2lu ", *u);
            }
        }
        fprintf(stderr, "\n");
    }

    // unlock coverage map
    if (conf->fast_mutex)
        pthread_mutex_unlock(&mutex_coverage);
}

void Alignment::alter_coverage_vector(vector<vector<vector<unsigned long> > > &cov_change, vector<vector<set<unsigned long> > > &genome_pos, bool is_curr_best) {

    // alter coverage iteratively in all coverage vectors
    for (size_t a = 0; a < cov_change.size(); a++) {
        // iterate over all segments of each coverage vector (e.g., caused by cut introns)
        for (size_t b = 0; b < cov_change.at(a).size(); b++) {

            if (conf->debug) {
                this->print();
                fprintf(stderr, "Cov before change:\n");
                for (vector<unsigned long>::iterator u = cov_change.at(a).at(b).begin(); u != cov_change.at(a).at(b).end(); u++) {
                    fprintf(stderr, "%2lu ", *u);
                }
                fprintf(stderr, "\n");
                for (set<unsigned long>::iterator v = genome_pos.at(a).at(b).begin(); v != genome_pos.at(a).at(b).end(); v++) {
                    fprintf(stderr, "%2lu ", *v);
                }
                fprintf(stderr, "\n");
            }

            // alter coverage where necessary
            set<unsigned long>::iterator gp;
            unsigned long chrm_pos = this->start;
            size_t curr_pos;

            for (size_t i = 0; i < this->sizes.size(); i++) {
                switch (this->operations.at(i)) {
                    case 'M': case 'D': {   for (int j = 0; j < this->sizes.at(i); j++) {
                                                gp = genome_pos.at(a).at(b).find(chrm_pos + genData->chr_size_cum.at(this->chr - 1));
                                                if (gp != genome_pos.at(a).at(b).end()) {
                                                    curr_pos = distance(genome_pos.at(a).at(b).begin(), gp);
                                                    if (this->is_best && is_curr_best) {
                                                        cov_change.at(a).at(b).at(curr_pos++) -= 1;
                                                    } else{
                                                        cov_change.at(a).at(b).at(curr_pos++) += 1;
                                                    }
                                                }
                                                chrm_pos++;
                                            }; 
                                            break;
                                        }
                    case 'N': { chrm_pos += this->sizes.at(i);
                                break;
                              }
                }
            }

            if (conf->debug) {
                fprintf(stderr, "Cov after change:\n");
                for (vector<unsigned long>::iterator u = cov_change.at(a).at(b).begin(); u != cov_change.at(a).at(b).end(); u++) {
                    fprintf(stderr, "%2li ", *u);
                }
                fprintf(stderr, "\n");
            }
        }
    }
}


void Alignment::update_coverage_map(bool positive) {

    vector<unsigned int>::iterator idx = genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].begin() + this->start;
    unsigned long pos = this->start;

    for (size_t i = 0; i < this->sizes.size(); i++) {
        switch (this->operations.at(i)) {
            case 'M': case 'D': for (int j = 0; j < this->sizes.at(i); j++) {
                                    if (idx < genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end()) {
                                        if (conf->fast_mutex)
                                            pthread_mutex_lock(&mutex_coverage);
                                        *idx += (*idx > 0 || positive) ? (2*positive - 1) : 0; 
                                        if (conf->fast_mutex)
                                            pthread_mutex_unlock(&mutex_coverage);
                                        idx++;
                                        pos++;
                                    }
                                }; 
                                break;
            case 'N': { if (genData->intron_coverage_map.find(pair<unsigned char, unsigned char>(this->chr, this->strand)) == genData->intron_coverage_map.end()) {
                            map<pair<unsigned long, unsigned long>, unsigned int> tmp;
                            genData->intron_coverage_map.insert(pair<pair<unsigned char, unsigned char>, map< pair<unsigned long, unsigned long>, unsigned int> >(pair<unsigned char, unsigned char>(this->chr, this->strand), tmp));
                        }
                        pthread_mutex_lock(&mutex_int_coverage);
                        map< pair<unsigned long, unsigned long>, unsigned int>::iterator it = genData->intron_coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].find(pair<unsigned long, unsigned long>(pos, pos + this->sizes.at(i) - 1));
                        if (it != genData->intron_coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end()) {
                            assert(it->second > 0 || positive);
                            it->second += (it->second > 0 || positive) ? (2*positive - 1) : 0;
                        } else
                            genData->intron_coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].insert(pair<pair<unsigned long, unsigned long>, unsigned int>(pair<unsigned long, unsigned long>(pos, pos + this->sizes.at(i) - 1), 1));
                        pthread_mutex_unlock(&mutex_int_coverage);
                        pos += this->sizes.at(i);
                        idx += this->sizes.at(i);
                      }
        }
        if (idx >= genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end())
            break;
    }
}

unsigned long Alignment::get_end() {
    // return alignment end in CLOSED intervals
    
    unsigned long end = this->start;
    for (size_t i = 0; i < this->sizes.size(); i++) {
        if (this->operations.at(i) == 'N' || this->operations.at(i) == 'D' || this->operations.at(i) == 'M') {
            end += this->sizes.at(i);
        }
    }
    
    return end == this->start?end:end - 1;
}

bool Alignment::compare_edit_ops(const Alignment &left, const Alignment &right) {
    return left.edit_ops < right.edit_ops;
}

bool Alignment::is_spliced() {
    return find(this->operations.begin(), this->operations.end(), 'N') != this->operations.end();
}

vector< pair<unsigned long, unsigned int> > Alignment::get_intron_coords() {
    
    vector< pair<unsigned long, unsigned int> > coords;
    unsigned long intron_start = this->start;

    for (size_t i = 0; i < this->operations.size(); i++) {
        if (this->operations.at(i) == 'M' || this->operations.at(i) == 'D') {
            intron_start += sizes.at(i);
        } else if (this->operations.at(i) == 'N') {
            coords.push_back(make_pair(intron_start, this->sizes.at(i)));
        }
    }
    return coords;
}

void Alignment::get_blocks(vector<pair<unsigned long, unsigned long> > &blocks) {
    
    unsigned long curr_pos = this->start;
    unsigned long last_pos = this->start;

    for (size_t i = 0; i < this->operations.size(); i++) {
        if (this->operations.at(i) == 'M' || this->operations.at(i) == 'D') {
            curr_pos += sizes.at(i);
        } else if (this->operations.at(i) == 'N') {
            blocks.push_back(make_pair(last_pos, curr_pos - 1));
            curr_pos += sizes.at(i);
            last_pos = curr_pos;
        }
    }
    blocks.push_back(make_pair(last_pos, curr_pos - 1));
}

set<unsigned long> Alignment::get_overlap(Alignment &other) {
    set<unsigned long> local_pos;
    set<unsigned long> overlap;

    if (this->chr != other.chr || this->strand != other.strand)
        return overlap;

    unsigned long genome_pos = this->start;
    for (size_t i = 0; i < this->operations.size(); i++) {
        if (this->operations.at(i) == 'M' || this->operations.at(i) == 'D') {
            for  (int j = 0; j < this->sizes.at(i); j++) {
                local_pos.insert(genome_pos++);
            }
        } else if (this->operations.at(i) == 'N') {
            genome_pos += this->sizes.at(i);
        }
    }

    genome_pos = other.start;
    for (size_t i = 0; i < other.operations.size(); i++) {
        if (other.operations.at(i) == 'M' || other.operations.at(i) == 'D') {
            for  (int j = 0; j < other.sizes.at(i); j++) {
                if (local_pos.find(genome_pos) != local_pos.end())
                    overlap.insert(genome_pos++);
            }
        } else if (other.operations.at(i) == 'N') {
            genome_pos += other.sizes.at(i);
        }
    }
    return overlap;
}


set<unsigned long> Alignment::get_exon_pos(unsigned int window_size) {

    set<unsigned long> position_set;

    unsigned long genome_pos = (window_size >= this->start) ? 0 : this->start - window_size;
    unsigned int offset = (conf->window_size >= this->start )? this->start : conf->window_size;
    
    // add global chromosome offset
    genome_pos += genData->chr_size_cum.at(this->chr - 1);

    // positions of preceding window
    for (unsigned int k = 0; k < offset; k++) {
        position_set.insert(genome_pos++);
    }

    // positions of actual alignment segments (exonic)
    for (size_t i = 0; i < this->operations.size(); i++) {
        if (this->operations.at(i) == 'M' || this->operations.at(i) == 'D') {
            for  (int j = 0; j < this->sizes.at(i); j++) {
                position_set.insert(genome_pos++);
            }
        } else if (this->operations.at(i) == 'N') {
            genome_pos += this->sizes.at(i);
        }
    }

    // positions of succeding window
    for (unsigned int k = 0; k < window_size; k++) {
        if ((genome_pos - genData->chr_size_cum.at(this->chr - 1)) < genData->coverage_map[ make_pair(this->chr, this->strand) ].size()) {
            position_set.insert(genome_pos++);
        } else {
            break;
        }
    }

    return position_set;
}


vector<set<unsigned long> > Alignment::get_genome_pos(unsigned int window_size) {

    vector<set<unsigned long> > position_set_vec;
    set<unsigned long> position_set;

    unsigned long genome_pos = (window_size >= this->start) ? 0 : this->start - window_size;
    unsigned int offset = (conf->window_size >= this->start )? this->start : conf->window_size;
    
    // add global chromosome offset
    unsigned long pos_offset = genData->chr_size_cum.at(this->chr - 1);
    genome_pos += pos_offset;

    pair<unsigned int, unsigned char> chr_strand = make_pair(this->chr, this->strand);

    // positions of preceding window
    for (unsigned int k = 0; k < offset; k++) {
        position_set.insert(genome_pos++);
        if (conf->use_brkpts && genome_pos - pos_offset < genData->breakpoint_map[ chr_strand ].size() && genData->breakpoint_map[ chr_strand ].at(genome_pos - pos_offset)) {
            position_set.clear();
        }
    }

    // positions of genomic alignment segments
    for (size_t i = 0; i < this->operations.size(); i++) {
        if (this->operations.at(i) == 'M' || this->operations.at(i) == 'D' || (this->operations.at(i) == 'N'  && this->sizes.at(i) <= (2 * (int) conf->window_size))) {
            for  (int j = 0; j < this->sizes.at(i); j++) {
                position_set.insert(genome_pos++);
                if (conf->use_brkpts && genome_pos - pos_offset < genData->breakpoint_map[ chr_strand ].size() && genData->breakpoint_map[ chr_strand ].at(genome_pos - pos_offset)) {
                    position_set_vec.push_back(position_set);
                    position_set.clear();
                }
            }
        } else if (this->operations.at(i) == 'N') {
            size_t pos2brkpt = 0;
            for  (size_t j = 0; j < conf->window_size; j++) {
                position_set.insert(genome_pos++);
                pos2brkpt++;
                if (conf->use_brkpts && genome_pos - pos_offset < genData->breakpoint_map[ chr_strand ].size() && genData->breakpoint_map[ chr_strand ].at(genome_pos - pos_offset)) {
                    position_set_vec.push_back(position_set);
                    position_set.clear();
                    break;
                }
            }
            genome_pos += (this->sizes.at(i) - conf->window_size - pos2brkpt);
            if (position_set.size() > 0) {
                position_set_vec.push_back(position_set);
                position_set.clear();
            }
            for  (size_t j = 0; j < conf->window_size; j++) {
                position_set.insert(genome_pos++);
                if (conf->use_brkpts && genome_pos - pos_offset < genData->breakpoint_map[ chr_strand ].size() && genData->breakpoint_map[ chr_strand ].at(genome_pos - pos_offset)) {
                    position_set.clear();
                }
            }
        }
    }

    // positions of succeding window
    for (unsigned int k = 0; k < window_size; k++) {
        if ((genome_pos - pos_offset) < genData->coverage_map[ chr_strand ].size()) {
            position_set.insert(genome_pos++);
            if (conf->use_brkpts && (genome_pos - pos_offset) < genData->breakpoint_map[ chr_strand ].size() && genData->breakpoint_map[ chr_strand ].at(genome_pos - pos_offset)) {
                position_set_vec.push_back(position_set);
                position_set.clear();
                break;
            }
        } else {
            break;
        }
    }

    if (position_set.size() > 0)
        position_set_vec.push_back(position_set);

    return position_set_vec;
}

void Alignment::clear() {
    this->chr = 0;
    this->start = 0;
    this->operations.clear();
    this->sizes.clear();
    this->is_best = false;
    this->edit_ops = 0;
    this->quality = 0;
    this->reversed = false;
    this->is_secondary = false;
    this->strand = '+';
}

void Alignment::print() {
    fprintf(stdout, "chr: %i, strand: %c start: %lu, cigar: ", this->chr, this->strand, this->start);
    for (unsigned int i = 0; i < this->sizes.size(); i++) {
        fprintf(stdout, "%i%c", this->sizes.at(i), this->operations.at(i));
    }
    fprintf(stdout, ", is_best: %i edit ops: %i\n", this->is_best, this->edit_ops);
}

void Alignment::determine_gaps() {
    unsigned int pos = this->start;

    for (size_t i = 0; i < this->operations.size(); i++) {
        if (this->operations.at(i) == 'M' || this->operations.at(i) == 'N')
            pos += sizes.at(i);
        else if (this->operations.at(i) == 'I') 
            this->insertions.insert(make_pair(pos, this->sizes.at(i)));
        else if (this->operations.at(i) == 'D') {
            for (int j = 0; j < this->sizes.at(i); j++)
                this->deletions.insert(pos + j);
            pos += this->sizes.at(i);
        }
    }
}

unsigned int Alignment::get_length() {
    unsigned int length = 0;
    for (size_t i = 0; i < this->operations.size(); i++) {
        if (this->operations.at(i) == 'M' || this->operations.at(i) == 'D') {
            length += this->sizes.at(i);
        }
    }
    return length;
}
