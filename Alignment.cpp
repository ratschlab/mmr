#include <algorithm>

#include "assert.h"
#include "Alignment.h"
#include "Utils.h"

extern GeneralData* genData;
extern Config* conf;

extern pthread_mutex_t mutex_coverage;
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
                id = id.substr(0, id.size() - 2);
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

void Alignment::fill_coverage_vector(vector<unsigned long> &cov_keep) {

    vector<unsigned int>::iterator cov_idx = genData->coverage_map[ make_pair(this->chr, this->strand) ].begin();
    size_t chrm_pos = 0;
    if (conf->window_size < this->start) {
        advance(cov_idx, this->start - conf->window_size);
        chrm_pos = (this->start - conf->window_size);
    }
//    size_t curr_pos = distance(genome_pos.begin(), genome_pos.find(chrm_pos + genData->chr_size_cum.at(this->chr - 1)));
    unsigned int offset = (conf->window_size >= this->start)?this->start:conf->window_size;

    // lock coverage map
    pthread_mutex_lock(&mutex_coverage);

    if (conf->debug) {
        this->print();
        fprintf(stderr, "cov keep before fill:\n");
        for (vector<unsigned long>::iterator u = cov_keep.begin(); u != cov_keep.end(); u++) {
            fprintf(stderr, "%2lu ", *u);
        }
        fprintf(stderr, "\n");
    }

    // get coverage from preceding windows
    for (size_t i = 0; i < offset; i++) {
        if (cov_idx < genData->coverage_map[ make_pair(this->chr, this->strand) ].end()) {
            cov_keep.push_back(*cov_idx);
            cov_idx++;
            chrm_pos++;
        } else {
            break;
        }
    }
    pthread_mutex_unlock(&mutex_coverage);

    // get coverage from alignment exons
    for (size_t i = 0; i < this->sizes.size(); i++) {
        switch (this->operations.at(i)) {
            case 'M': case 'D': {   //curr_pos = distance(genome_pos.begin(), genome_pos.find(chrm_pos + genData->chr_size_cum.at(this->chr - 1)));
                                    pthread_mutex_lock(&mutex_coverage);
                                    for (int j = 0; j < this->sizes.at(i); j++) {
                                        if (cov_idx < genData->coverage_map[ make_pair(this->chr, this->strand) ].end()) {
                                            cov_keep.push_back(*cov_idx);
                                            cov_idx++;
                                            chrm_pos++;
                                        }
                                    }; 
                                    pthread_mutex_unlock(&mutex_coverage);
                                    break;
                                }
            case 'N': { cov_idx += this->sizes.at(i);
                        chrm_pos += this->sizes.at(i);
                        break;
                      }
        }
    }

    // get coverage from following windows
    pthread_mutex_lock(&mutex_coverage);
    for (size_t i = 0; i < conf->window_size; i++) {
        if (cov_idx < genData->coverage_map[ make_pair(this->chr, this->strand) ].end()) {
            cov_keep.push_back(*cov_idx);
            cov_idx++;
        }
        else {
            break;
        }
    }

    if (conf->debug) {
        fprintf(stderr, "cov keep after fill:\n");
        for (vector<unsigned long>::iterator u = cov_keep.begin(); u != cov_keep.end(); u++) {
            fprintf(stderr, "%2lu ", *u);
        }
        fprintf(stderr, "\n");
    }

    // unlock coverage map
    pthread_mutex_unlock(&mutex_coverage);
}

void Alignment::alter_coverage_vector(vector<vector<unsigned long> > &cov_change, vector<set<unsigned long> > &genome_pos, bool is_curr_best) {

    // alter coverage iteratively in all coverage vectors
    for (size_t a = 0; a < cov_change.size(); a++) {

        if (conf->debug) {
            this->print();
            fprintf(stderr, "Cov before change:\n");
            for (vector<unsigned long>::iterator u = cov_change.at(a).begin(); u != cov_change.at(a).end(); u++) {
                fprintf(stderr, "%2lu ", *u);
            }
            fprintf(stderr, "\n");
            for (set<unsigned long>::iterator v = genome_pos.at(a).begin(); v != genome_pos.at(a).end(); v++) {
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
                                            gp = genome_pos.at(a).find(chrm_pos + genData->chr_size_cum.at(this->chr - 1));
                                            if (gp != genome_pos.at(a).end()) {
                                                curr_pos = distance(genome_pos.at(a).begin(), gp);
                                                if (this->is_best && is_curr_best) {
                                                    cov_change.at(a).at(curr_pos++) -= 1;
                                                } else{
                                                    cov_change.at(a).at(curr_pos++) += 1;
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
            for (vector<unsigned long>::iterator u = cov_change.at(a).begin(); u != cov_change.at(a).end(); u++) {
                fprintf(stderr, "%2li ", *u);
            }
            fprintf(stderr, "\n");
        }
    }
}

//void Alignment::fill_coverage_vectors(vector<unsigned long> &cov_keep, vector<unsigned long> &cov_change, vector<bool> &intronic, unsigned long first_start, bool is_curr_best) {
void Alignment::fill_coverage_vectors(vector<unsigned long> &cov_keep, vector<unsigned long> &cov_change, set<unsigned int> &genome_pos, unsigned long first_start, bool is_curr_best) {

    // cov_keep -> coverage, if we keep the current best assignment
    // cov_change -> coverage, if we change the current best assignment

    vector<unsigned int>::iterator cov_idx = genData->coverage_map[ make_pair(this->chr, this->strand) ].begin();

    unsigned int offset = (conf->window_size >= this->start)?this->start:conf->window_size;
    // position in global coverage coordinates (for all chromosomes)
    size_t curr_pos = distance(genome_pos.begin(), genome_pos.find(this->start + genData->chr_size_cum.at(this->chr - 1)));

    // set start coordinates
    if (conf->window_size < this->start) { 
        //fprintf(stdout, "this start: %lu\n", this->start);
        //fprintf(stdout, "window %i\n", conf->window_size);
        //fprintf(stdout, "diff %lu\n", this->start - conf->window_size);
        advance(cov_idx, this->start - conf->window_size);
        //curr_pos += (this->start - conf->window_size - first_start);
        //fprintf(stdout, "first start %lu\n", first_start);
        //fprintf(stdout, "curr_pos %lu\n", curr_pos);
    }

    // lock coverage map
    pthread_mutex_lock(&mutex_coverage);
    //fprintf(stdout, "Coverage:\n");
    //for (vector<unsigned int>::iterator tt = genData->coverage_map[ make_pair(this->chr, this->strand) ].begin(); tt != genData->coverage_map[ make_pair(this->chr, this->strand) ].end(); tt++) {
    //    fprintf(stdout, "%u ", *tt);
    //}
    //fprintf(stdout, "\n\n");

    // get coverage from preceding windows
    for (size_t i = 0; i < offset; i++) {
        if (cov_idx < genData->coverage_map[ make_pair(this->chr, this->strand) ].end()) {

            assert(cov_keep.size() == cov_change.size());

            if (curr_pos >= cov_keep.size()) {
                //fprintf(stdout, "curr pos: %i\n", curr_pos);
                cov_keep.push_back(*cov_idx);
                cov_change.push_back(*cov_idx);
            }
            cov_idx++;
            curr_pos++;
        } else {
            break;
        }
    }

    // get coverage from alignments
    for (size_t i = 0; i < this->sizes.size(); i++) {
        switch (this->operations.at(i)) {
            case 'M': case 'D': { for (int j = 0; j < this->sizes.at(i); j++) {
                                    if (cov_idx < genData->coverage_map[ make_pair(this->chr, this->strand) ].end()) {
                                        if (this->is_best && is_curr_best) {
                                            if (curr_pos < cov_change.size()) {
                                                //fprintf(stdout, "changed %lu to %lu\n", cov_change.at(curr_pos), cov_change.at(curr_pos) - 1);
                                                cov_change.at(curr_pos) -= 1;
                                            } else {
                                                assert(*cov_idx > 0);
                                                cov_change.push_back((*cov_idx) - 1);
                                                cov_keep.push_back(*cov_idx);
                                            }
                                        } else{
                                            if (curr_pos < cov_change.size()) {
                                                cov_change.at(curr_pos) += 1;
                                            } else {
                                                cov_change.push_back((*cov_idx) + 1);
                                                cov_keep.push_back(*cov_idx);
                                            }
                                        }
                                        //fprintf(stdout, "curr pos: %i\n", curr_pos);
                                        cov_idx++;
                                        curr_pos++;
                                     }
                                  }; 
                                  break;
                                }
            case 'N': { cov_idx += this->sizes.at(i);
                        break;
                      }
        }
    }

    // get coverage from following windows
    for (size_t i = 0; i < conf->window_size; i++) {
        if (cov_idx < genData->coverage_map[ make_pair(this->chr, this->strand) ].end()) {

            assert(cov_keep.size() == cov_change.size());

            if (curr_pos >= cov_change.size()) {
                cov_change.push_back(*cov_idx);
                cov_keep.push_back(*cov_idx);
                //fprintf(stdout, "curr pos: %i\n", curr_pos);
            }
            cov_idx++;
            curr_pos++;
        }
        else {
            break;
        }
    }
    // unlock coverage map
    pthread_mutex_unlock(&mutex_coverage);

    /*fprintf(stdout, "keep:\n");
    for (size_t uu = 0; uu < cov_keep.size(); uu++) {
        fprintf(stdout, "%lu ", cov_keep.at(uu));
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "change:\n");
    for (size_t uu = 0; uu < cov_change.size(); uu++) {
        fprintf(stdout, "%lu ", cov_change.at(uu));
    }
    fprintf(stdout, "\n\n");
    */
}

void Alignment::update_coverage_map(bool positive) {

    //pthread_mutex_lock(&mutex_coverage);
    vector<unsigned int>::iterator idx = genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].begin() + this->start;
    unsigned long pos = this->start;

    for (size_t i = 0; i < this->sizes.size(); i++) {
        switch (this->operations.at(i)) {
            case 'M': case 'D': for (int j = 0; j < this->sizes.at(i); j++) {
                                    if (idx < genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end()) {
                                        pthread_mutex_lock(&mutex_coverage);
                                        *idx += (*idx > 0 || positive) ? (2*positive - 1) : 0; 
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
                        map< pair<unsigned long, unsigned long>, unsigned int>::iterator it = genData->intron_coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].find(pair<unsigned long, unsigned long>(pos, pos + this->sizes.at(i) - 1));
                        if (it != genData->intron_coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end()) {
                            assert(it->second > 0 || positive);
                            it->second += (it->second > 0 || positive) ? (2*positive - 1) : 0;
                        } else
                            genData->intron_coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].insert(pair<pair<unsigned long, unsigned long>, unsigned int>(pair<unsigned long, unsigned long>(pos, pos + this->sizes.at(i) - 1), 1));
                        pos += this->sizes.at(i);
                        idx += this->sizes.at(i);
                      }
        }
        if (idx >= genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end())
            break;
    }

    //pthread_mutex_unlock(&mutex_coverage);
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
/*    bool spliced = false;
    for (vector<char>::iterator it = this->operations.begin(); it != this->operations.end(); it++) {
        if ((*it) == 'N') {
            spliced = true;
            break;
        }
    }*/
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

set<unsigned long> Alignment::get_genome_pos(unsigned int window_size) {

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
