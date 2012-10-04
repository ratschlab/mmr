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

string Alignment::fill(char* sl, unsigned char &pair) {

    int idx = 0;
    string id;

    while (sl != NULL) {
        if (idx == 0) { 
            id = sl;
        } else if (idx == 1) {
            pair = (atoi(sl) & 128);
            this->reversed = ((atoi(sl) & 16) == 16);
        } else if (idx == 2) {
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
            break;
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

pair<double, double> Alignment::get_variance_loss(set<unsigned long> covered_pos, set<unsigned long> not_covered_pos, set<unsigned long> mate_covered_pos, bool mate_is_best) {

    vector<unsigned int>::iterator cov_idx = genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].begin();
    vector<unsigned int>::iterator intron_cov_idx;
    
    vector<unsigned int> exon_cov_with;
    vector<unsigned int> exon_cov_without;
    vector<unsigned int> intron_cov;

    unsigned long genome_pos = 0;
    unsigned int adjust_up = 0;
    unsigned int adjust_down = 0;
    unsigned int adjust_mate = 0;
    unsigned int offset = (conf->window_size >= this->start)?this->start:conf->window_size;
    size_t step_size = 1;

    // set start coordinates
    if (conf->window_size < this->start) {  
        cov_idx += (this->start - conf->window_size);
        genome_pos += (this->start - conf->window_size);
    }
        
    // get coverage from preceding windows
    for (size_t i = 0; i < offset; i++) {
        if (cov_idx < genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end()) {
            adjust_up = (not_covered_pos.find(genome_pos) != not_covered_pos.end()) ? 1 : 0;
            adjust_down = (covered_pos.find(genome_pos) != covered_pos.end()) ? 1 : 0;
            adjust_mate = (mate_covered_pos.find(genome_pos) != mate_covered_pos.end()) ? 1 : 0;
            if (this->is_best && mate_is_best) {
                exon_cov_with.push_back((*cov_idx) - adjust_down);
                exon_cov_without.push_back((*cov_idx) + adjust_up - adjust_mate);
            } else {
                adjust_mate = mate_is_best ? 0 : adjust_mate;
                exon_cov_with.push_back((*cov_idx) - adjust_down + adjust_mate);
                exon_cov_without.push_back((*cov_idx) + adjust_up);
            }
            cov_idx++;
            genome_pos++;
        } else {
            break;
        }
    }

    // get coverage from alignments
    for (size_t i = 0; i < this->sizes.size(); i++) {
        switch (this->operations.at(i)) {
            case 'M': case 'D': { for (int j = 0; j < this->sizes.at(i); j++) {
                                    if (cov_idx < genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end()) {
                                        adjust_up = (not_covered_pos.find(genome_pos) != not_covered_pos.end()) ? 1 : 0;
                                        adjust_down = (covered_pos.find(genome_pos) != covered_pos.end()) ? 1 : 0;

                                        if (this->is_best && mate_is_best) {
                                            exon_cov_with.push_back(*cov_idx);
                                            exon_cov_without.push_back(max((unsigned int) 0, (*cov_idx) - 1 + adjust_up));
                                        } else if (this->is_best) { // it means this is the candidate in the pair case and does not change
                                            exon_cov_with.push_back(*cov_idx);
                                            exon_cov_without.push_back(*cov_idx);
                                        } else {
                                            exon_cov_with.push_back((*cov_idx) + 1 - adjust_down);
                                            exon_cov_without.push_back(*cov_idx);
                                        }
                                        cov_idx++;
                                        genome_pos++;
                                     }
                                  }; 
                                  break;
                                }
            case 'N': { step_size = max(this->sizes.at(i) / 50, 1);
                        for (int j = conf->intron_offset; j < this->sizes.at(i) - conf->intron_offset; j += step_size) {
                            if (intron_cov_idx < genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end()) {
                                adjust_down = (covered_pos.find(genome_pos + j) != covered_pos.end()) ? 1 : 0;
                                intron_cov_idx = cov_idx + j;
                                intron_cov.push_back((*intron_cov_idx) - adjust_down); 
                            }
                        }; 
                        cov_idx += this->sizes.at(i);
                        genome_pos += this->sizes.at(i);
                        break;
                      }
        }
    }

    // get coverage from following windows
    for (size_t i = 0; i < conf->window_size; i++) {
        if (cov_idx < genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end()) {
            adjust_up = (not_covered_pos.find(genome_pos) != not_covered_pos.end()) ? 1 : 0;
            adjust_down = (covered_pos.find(genome_pos) != covered_pos.end()) ? 1 : 0;
            adjust_mate = (mate_covered_pos.find(genome_pos) != mate_covered_pos.end()) ? 1 : 0;
            if (this->is_best && mate_is_best) {
                exon_cov_with.push_back((*cov_idx) - adjust_down);
                exon_cov_without.push_back((*cov_idx) + adjust_up - adjust_mate);
            } else {
                adjust_mate = mate_is_best ? 0 : adjust_mate;
                exon_cov_with.push_back((*cov_idx) - adjust_down + adjust_mate);
                exon_cov_without.push_back((*cov_idx) + adjust_up);
            }
            cov_idx++;
            genome_pos++;
        }
        else {
            break;
        }
    }
  
    // compute variance loss
    vector<unsigned int> empty_cov;
    double loss_with = get_variance(exon_cov_with, intron_cov);
    double loss_without = get_variance(exon_cov_without, empty_cov);

    return pair<double, double>(loss_with, loss_without);
}

void Alignment::update_coverage_map(bool positive) {

    pthread_mutex_lock(&mutex_coverage);
    //fprintf(stdout, "update coverage map on chr %i\n", this->chr);
    vector<unsigned int>::iterator idx = genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].begin() + this->start;
    unsigned long pos = this->start;

    for (size_t i = 0; i < this->sizes.size(); i++) {
        switch (this->operations.at(i)) {
            case 'M': case 'D': for (int j = 0; j < this->sizes.at(i); j++) {
                                    if (idx < genData->coverage_map[pair<unsigned char, unsigned char>(this->chr, this->strand)].end()) {
                                        *idx += (*idx > 0 || positive) ? (2*positive - 1) : 0; 
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

    pthread_mutex_unlock(&mutex_coverage);
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

bool Alignment::comparator(const Alignment &left, const Alignment &right) {
    return left.edit_ops < right.edit_ops;
}

bool Alignment::is_spliced() {
    bool spliced = false;
    for (vector<char>::iterator it = this->operations.begin(); it != this->operations.end(); it++) {
        if ((*it) == 'N') {
            spliced = true;
            break;
        }
    }
    return spliced;
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

    /*if (window_size > 0) {
        for (unsigned long i = max(0, this->start - window_size); i < this->start; i++) {
            position_set.insert(i);
        }
    }*/
    unsigned long genome_pos = this->start;
    for (size_t i = 0; i < this->operations.size(); i++) {
        if (this->operations.at(i) == 'M' || this->operations.at(i) == 'D') {
            for  (int j = 0; j < this->sizes.at(i); j++) {
                position_set.insert(genome_pos++);
            }
        } else if (this->operations.at(i) == 'N') {
            genome_pos += this->sizes.at(i);
        }
    }
    /*if (window_size > 0) {
        for (unsigned long i = genome_pos; i < genData->chr_size.at(this->chr); i++) {
            position_set.insert(i);
        }
    }*/

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
    this->strand = '+';
}

void Alignment::print() {
    fprintf(stdout, "chr: %i, strand: %c start: %lu, cigar: ", this->chr, this->strand, this->start);
    for (unsigned int i = 0; i < this->sizes.size(); i++) {
        fprintf(stdout, "%i%c", this->sizes.at(i), this->operations.at(i));
    }
    fprintf(stdout, ", is_best: %i\n", this->is_best);
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
