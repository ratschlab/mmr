
#include <assert.h>
#include <limits>

#include "OnlineData.h"
#include "Utils.h"
#include "config.h"

extern Config* conf;

extern pthread_mutex_t mutex_coverage;
extern pthread_mutex_t mutex_fifo;
extern pthread_mutex_t mutex_done;
extern pthread_mutex_t mutex_best_left;
extern pthread_mutex_t mutex_best_right;
extern pthread_mutex_t mutex_counter;

void OnlineData::process_data_online(GeneralData* genData) {

    vector<vector<Alignment>::iterator> active_left_reads;
    vector<vector<Alignment>::iterator> active_right_reads;

    vector<vector<Alignment>::iterator>::iterator lv_idx;
    vector<vector<Alignment>::iterator>::iterator rv_idx;
    vector<vector<Alignment>::iterator>::iterator curr_best;

    vector<Alignment>::iterator best_left_idx;
    vector<Alignment>::iterator best_right_idx;

    set<vector<Alignment>::iterator> ignore_idx_left;
    set<vector<Alignment>::iterator> ignore_idx_right;

    unsigned int num_changed = 0;
    double min_loss = numeric_limits<double>::max();
    double total_min_loss = 0.0;
    double loss = 0.0;
    bool found_pairs = false;

    if (conf->pre_filter) {
        ignore_idx_left = filter_alignments(this->left_reads);
        pthread_mutex_lock(&mutex_best_left);
        if (ignore_idx_left.find(this->left_reads.begin() + genData->best_left[this->last_id]) != ignore_idx_left.end()) {
            vector<Alignment>::iterator curr_align = (this->left_reads.begin() + genData->best_left[this->last_id]);
            curr_align->is_best = false;
            curr_align->update_coverage_map(0);
            for (vector<Alignment>::iterator v_idx = this->left_reads.begin(); v_idx != this->left_reads.end(); v_idx++) {
                if (ignore_idx_left.find(v_idx) == ignore_idx_left.end()) {
                    genData->best_left[this->last_id] = v_idx - this->left_reads.begin();
                    v_idx->is_best = true;
                    v_idx->update_coverage_map(1);
                    break;
                }
            }
        }
        pthread_mutex_unlock(&mutex_best_left);

        ignore_idx_right = filter_alignments(this->right_reads);
        pthread_mutex_lock(&mutex_best_right);
        if (ignore_idx_right.find(this->right_reads.begin() + genData->best_right[this->last_id]) != ignore_idx_right.end()) {
            vector<Alignment>::iterator curr_align = (this->right_reads.begin() + genData->best_right[this->last_id]);
            curr_align->is_best = false;
            curr_align->update_coverage_map(0);
            for (vector<Alignment>::iterator v_idx = this->right_reads.begin(); v_idx != this->right_reads.end(); v_idx++) {
                if (ignore_idx_right.find(v_idx) == ignore_idx_right.end()) {
                    genData->best_right[this->last_id] = v_idx - this->right_reads.begin();
                    v_idx->is_best = true;
                    v_idx->update_coverage_map(1);
                    break;
                }
            }
        }
        pthread_mutex_unlock(&mutex_best_right);
    }

    get_active_reads(this->last_id, ignore_idx_left, ignore_idx_right, active_left_reads, active_right_reads, genData, found_pairs);

    if (found_pairs) {
        bool best_found = false;
        for(lv_idx = active_left_reads.begin(), rv_idx = active_right_reads.begin(); lv_idx < active_left_reads.end() && rv_idx  < active_right_reads.end(); lv_idx++, rv_idx++) {
            if ((*rv_idx)->is_best && (*lv_idx)->is_best) {
                best_left_idx = (*lv_idx);
                best_right_idx = (*rv_idx);
                best_found = true;
                break;
            }
        }
        assert(best_found);

        bool changed = false;
        //double new_total_loss = 0.0;
        for(lv_idx = active_left_reads.begin(), rv_idx = active_right_reads.begin(); lv_idx != active_left_reads.end() && rv_idx != active_right_reads.end(); lv_idx++, rv_idx++) {
            if ((*lv_idx) == best_left_idx && (*rv_idx) == best_right_idx)
                continue;
            if (compare_pair(*lv_idx, *rv_idx, best_left_idx, best_right_idx, loss)) {

                changed = true;

                vector<Alignment>::iterator old_best_left = best_left_idx;
                vector<Alignment>::iterator old_best_right = best_right_idx;

                best_left_idx->is_best = false;
                best_right_idx->is_best = false;
                best_left_idx->update_coverage_map(0);
                best_right_idx->update_coverage_map(0);
                best_left_idx = *lv_idx;
                best_right_idx = *rv_idx;
                best_left_idx->is_best = true;
                best_right_idx->is_best = true;
                best_left_idx->update_coverage_map(1);
                best_right_idx->update_coverage_map(1);
                
                /*if (conf->iteration > 0) {
                    new_total_loss = genData->segments.get_total_loss();
                    if (new_total_loss > conf->last_loss) {
                        // do rollback
                        best_left_idx->is_best = false;
                        best_right_idx->is_best = false;
                        best_left_idx->update_coverage_map(0);
                        best_right_idx->update_coverage_map(0);
                        best_left_idx = old_best_left;
                        best_right_idx = old_best_right;
                        best_left_idx->is_best = true;
                        best_right_idx->is_best = true;
                        best_left_idx->update_coverage_map(1);
                        best_right_idx->update_coverage_map(1);
                        compare_pair(*lv_idx, *rv_idx, best_left_idx, best_right_idx, loss, true);
                        fprintf(stdout, "loss so far: %f\n", conf->last_loss);
                        fprintf(stdout, "loss now: %f\n\n", new_total_loss);
                        // roll-roll-back
                        best_left_idx->is_best = false;
                        best_right_idx->is_best = false;
                        best_left_idx->update_coverage_map(0);
                        best_right_idx->update_coverage_map(0);
                        best_left_idx = *lv_idx;
                        best_right_idx = *rv_idx;
                        best_left_idx->is_best = true;
                        best_right_idx->is_best = true;
                        best_left_idx->update_coverage_map(1);
                        best_right_idx->update_coverage_map(1);
                    }
                    conf->last_loss = new_total_loss;
                }*/

                pthread_mutex_lock(&mutex_best_left);
                genData->best_left[this->last_id] = (best_left_idx - this->left_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_left);
                pthread_mutex_lock(&mutex_best_right);
                genData->best_right[this->last_id] = (best_right_idx - this->right_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_right);
            }
            if (loss >= 0.0)
                min_loss = min(min_loss, loss);
        }
        total_min_loss += (min_loss < numeric_limits<double>::max()) ? min_loss : 0;
        if (changed) 
            num_changed++;

    } else {
        bool best_found = false;
        for (lv_idx = active_left_reads.begin(); lv_idx != active_left_reads.end(); lv_idx++) {
            if ((*lv_idx)->is_best) {
                curr_best = lv_idx;
                best_found = true;
                break;
            }
        }
        if (active_left_reads.size() > 0) {
            if (!best_found)
                fprintf(stderr, "The given input file is most probably not sorted by read-ID! Bailing out. \n");
            assert(best_found);
        }
        bool changed = false;
        for (lv_idx = active_left_reads.begin(); lv_idx != active_left_reads.end(); lv_idx++) {
            if (*lv_idx == *curr_best)
                continue;

            vector<Alignment>::iterator old_best = *curr_best;
            // check if lv_idx < curr_best
            if (compare_single(*lv_idx, *curr_best, loss)) {
                changed = true;
                (*curr_best)->is_best = false;
                (*curr_best)->update_coverage_map(0);
                *curr_best = *lv_idx;          
                (*curr_best)->is_best = true;
                (*curr_best)->update_coverage_map(1);
                /*if (conf->iteration > 0) {
                    double new_total_loss = genData->segments.get_total_loss();
                    if (new_total_loss > conf->last_loss) {
                        // do rollback
                        old_best->print();
                        (*curr_best)->print();
                        (*lv_idx)->print();
                        (*curr_best)->is_best = false;
                        (*curr_best)->update_coverage_map(0);
                        (*curr_best) = old_best;
                        (*curr_best)->is_best = true;
                        (*curr_best)->update_coverage_map(1);
                        (*lv_idx)->print();
                        (*curr_best)->print();
                        compare_single(*lv_idx, *curr_best, loss, true);
                        fprintf(stdout, "left single:\n");
                        fprintf(stdout, "loss so far: %f\n", conf->last_loss);
                        fprintf(stdout, "loss now: %f\n\n", new_total_loss);
                        // roll-roll-back
                        (*curr_best)->is_best = false;
                        (*curr_best)->update_coverage_map(0);
                        *curr_best = *lv_idx;          
                        (*curr_best)->is_best = true;
                        (*curr_best)->update_coverage_map(1);
                    }
                    conf->last_loss = new_total_loss;
                }*/

                pthread_mutex_lock(&mutex_best_left);
                genData->best_left[this->last_id] = (*curr_best - this->left_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_left);

            }
            if (loss >= 0.0)
                min_loss = min(min_loss, loss);
        }
        total_min_loss += (min_loss < numeric_limits<double>::max()) ? min_loss : 0;
        if (changed) 
            num_changed++;

        best_found = false;
        for (rv_idx = active_right_reads.begin(); rv_idx != active_right_reads.end(); rv_idx++) {
            if ((*rv_idx)->is_best) {
                curr_best = rv_idx;
                best_found = true;
                break;
            }
        }
        if (active_right_reads.size() > 0) {
            if (!best_found)
                fprintf(stderr, "The given input file is most probably not sorted by read-ID! Bailing out. \n");
            assert(best_found);
        }
        changed = false;
        min_loss = numeric_limits<double>::max();
        for (rv_idx = active_right_reads.begin(); rv_idx != active_right_reads.end(); rv_idx++) {
            if (*rv_idx == *curr_best)
                continue;
            
            vector<Alignment>::iterator old_best = *curr_best;
            // check if rv_idx < curr_best
            if (compare_single(*rv_idx, *curr_best, loss)) {
                changed = true;
                (*curr_best)->is_best = false;
                (*curr_best)->update_coverage_map(0);
                *curr_best = *rv_idx;          
                (*curr_best)->is_best = true;
                (*curr_best)->update_coverage_map(1);
                /*if (conf->iteration > 0) {
                    double new_total_loss = genData->segments.get_total_loss();
                    if (new_total_loss > conf->last_loss) {
                        old_best->print();
                        (*curr_best)->print();
                        // do rollback
                        (*curr_best)->is_best = false;
                        (*curr_best)->update_coverage_map(0);
                        (*curr_best) = old_best;
                        (*curr_best)->is_best = true;
                        (*curr_best)->update_coverage_map(1);
                        compare_single(*rv_idx, *curr_best, loss, true);
                        fprintf(stdout, "right single:\n");
                        fprintf(stdout, "loss so far: %f\n", conf->last_loss);
                        fprintf(stdout, "loss now: %f\n\n", new_total_loss);
                        // roll-roll-back
                        (*curr_best)->is_best = false;
                        (*curr_best)->update_coverage_map(0);
                        *curr_best = *rv_idx;          
                        (*curr_best)->is_best = true;
                        (*curr_best)->update_coverage_map(1);
                    }
                    conf->last_loss = new_total_loss;
                }*/

                pthread_mutex_lock(&mutex_best_right);
                genData->best_right[this->last_id] = (*curr_best - this->right_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_right);
            }
            if (loss >= 0.0)
                min_loss = min(min_loss, loss);
        }
        total_min_loss += (min_loss < numeric_limits<double>::max()) ? min_loss : 0;
        if (changed) 
            num_changed++;

    }
    pthread_mutex_lock(&mutex_counter);
    genData->total_loss += total_min_loss;
    genData->num_altered += num_changed;
    pthread_mutex_unlock(&mutex_counter);
    delete this;
}

void OnlineData::get_active_reads(string read_id, set<vector<Alignment>::iterator> &ignore_reads_left, set<vector<Alignment>::iterator> &ignore_reads_right, vector<vector<Alignment>::iterator> &active_left_reads, vector<vector<Alignment>::iterator> &active_right_reads, GeneralData* genData, bool &found_pairs) {

    // if pair processing, identify all compatible pairs
    // pairs are compatible, if they show same chr, opposite strands and
    // have an inner distance within the insert size range

    active_left_reads.clear();
    active_right_reads.clear();

    if (conf->use_pair_info && ! (this->right_reads.size() == 0 || this->left_reads.size() == 0)) {
        bool best_pair = false;
        for (vector<Alignment>::iterator lv_idx = this->left_reads.begin(); lv_idx != this->left_reads.end(); lv_idx++) {
            if (conf->pre_filter && (ignore_reads_left.find(lv_idx) != ignore_reads_left.end()))
                continue;
            for (vector<Alignment>::iterator rv_idx = this->right_reads.begin(); rv_idx != this->right_reads.end(); rv_idx++) {
                if (conf->pre_filter && (ignore_reads_right.find(rv_idx) != ignore_reads_right.end()))
                    continue;
                if (pair_is_valid(rv_idx, lv_idx)) {
                    best_pair = best_pair ? best_pair : (lv_idx->is_best && rv_idx->is_best);
                    active_left_reads.push_back(lv_idx);
                    active_right_reads.push_back(rv_idx);
                    found_pairs = true;
                }
            }
        }

        // no active pair is best alignment
        if (! best_pair && found_pairs) {
            for (vector<Alignment>::iterator lv_idx = this->left_reads.begin(); lv_idx != this->left_reads.end(); lv_idx++) {
                if (lv_idx->is_best) {
                    lv_idx->is_best = false;
                    lv_idx->update_coverage_map(0);
                    break;
                }
            }
            active_left_reads.front()->is_best = true;
            active_left_reads.front()->update_coverage_map(1);

            pthread_mutex_lock(&mutex_best_left);
            genData->best_left[read_id] = (active_left_reads.front() - this->left_reads.begin());
            pthread_mutex_unlock(&mutex_best_left);

            for (vector<Alignment>::iterator rv_idx = this->right_reads.begin(); rv_idx != this->right_reads.end(); rv_idx++) {
                if (rv_idx->is_best) {
                    rv_idx->is_best = false;
                    rv_idx->update_coverage_map(0);
                    break;
                }
            }
            active_right_reads.front()->is_best = true;
            active_right_reads.front()->update_coverage_map(1);

            pthread_mutex_lock(&mutex_best_right);
            assert(active_right_reads.front() - this->right_reads.begin() >= 0);
            genData->best_right[read_id] = (active_right_reads.front() - this->right_reads.begin());
            pthread_mutex_unlock(&mutex_best_right);
        }
        // did not find valid pairs
        if (! found_pairs) {
            active_left_reads.clear();
            for (vector<Alignment>::iterator lv_idx = this->left_reads.begin(); lv_idx != this->left_reads.end(); lv_idx++) {
                if (conf->pre_filter && (ignore_reads_left.find(lv_idx) != ignore_reads_left.end()))
                    continue;
                active_left_reads.push_back(lv_idx);
            }
            active_right_reads.clear();
            for (vector<Alignment>::iterator rv_idx = this->right_reads.begin(); rv_idx != this->right_reads.end(); rv_idx++) {
                if (conf->pre_filter && (ignore_reads_right.find(rv_idx) != ignore_reads_right.end()))
                    continue;
                active_right_reads.push_back(rv_idx);
            }
        }
    } else {
        for (vector<Alignment>::iterator lv_idx = this->left_reads.begin(); lv_idx != this->left_reads.end(); lv_idx++) {
            if (conf->pre_filter && (ignore_reads_left.find(lv_idx) != ignore_reads_left.end()))
                continue;
            active_left_reads.push_back(lv_idx);
        }
        for (vector<Alignment>::iterator rv_idx = this->right_reads.begin(); rv_idx != this->right_reads.end(); rv_idx++) {
            if (conf->pre_filter && (ignore_reads_right.find(rv_idx) != ignore_reads_right.end()))
                continue;
            active_right_reads.push_back(rv_idx);
        }
    }
}

char* OnlineData::parse_file(FILE* infile, char* last_line, GeneralData* genData, unsigned int &counter) {

    char line[1000] ;
    char cp_line[1000];
    char* ret = last_line;

    unsigned char pair_info = 0;

    Alignment curr_alignment;
    string id;
    this->last_id.clear();

    unordered_map <string, size_t, hash<string> >::iterator b_idx;

    this->left_reads.clear();
    this->right_reads.clear();

    while (true) {
        if (strlen(last_line) > 0 && strcmp(last_line, "samtools subprocess for reading terminated successfully\n")) { 
            strcpy(line, last_line);
            strcpy(last_line, "");
        }
        else {
            ret = fgets(line, sizeof(line), infile);
            counter++;

            if (conf->verbose && counter % 100000 == 0) 
                fprintf(stdout, "\n\t%i...", counter);
        }

        strcpy(cp_line, line);
        if (!ret) {
            break;
        }

        char* sl = strtok(line, "\t");

        curr_alignment.clear();
        id = curr_alignment.fill(sl, pair_info);

        if (id.size() == 0) {
            if (strcmp(cp_line, "samtools subprocess for reading terminated successfully\n")) {
                fprintf(stderr, "\nWARNING: SAM line incomplete! Ignoring line:\n%s\n", cp_line);
            } else {
                if (conf->verbose)
                    fprintf(stdout, "\n%s", cp_line);
            }
            continue ;
        }

        if (pair_info == 0) {
            // id == last_id or last_id is empty
            if (id.compare(this->last_id) == 0 || this->last_id.size() == 0) {
                b_idx = genData->best_left.find(id);
                if (b_idx == genData->best_left.end()) {
                    pthread_mutex_lock(&mutex_best_left);
                    genData->best_left.insert(pair<string, size_t>(id, this->left_reads.size()));
                    pthread_mutex_unlock(&mutex_best_left);
                    curr_alignment.is_best = true;
                    curr_alignment.update_coverage_map(1);
                }
                else if (b_idx->second == this->left_reads.size()) {
                    curr_alignment.is_best = true;
                    // TODO unsicher
                }
                this->left_reads.push_back(curr_alignment);
                this->last_id = id;
            } else {
                break ;
            }
        } else {
            // id == last_id or last_id is empty
            if ((! id.compare(this->last_id)) || this->last_id.size() == 0) {
                b_idx = genData->best_right.find(id);
                if (b_idx == genData->best_right.end()) {
                    pthread_mutex_lock(&mutex_best_right);
                    genData->best_right.insert(pair<string, size_t>(id, this->right_reads.size()));
                    pthread_mutex_unlock(&mutex_best_right);
                    curr_alignment.is_best = true;
                    curr_alignment.update_coverage_map(1);
                }
                else if (b_idx->second == this->right_reads.size()) {
                    curr_alignment.is_best = true;
                    // TODO unsicher
                }
                this->right_reads.push_back(curr_alignment);
                this->last_id = id;
            } else {
                break ;
            }
        }
    }

    strcpy(last_line, cp_line);
    return ret;
}

