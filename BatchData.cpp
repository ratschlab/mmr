#include "assert.h"

#include "BatchData.h"
#include "GeneralData.h"
#include "Alignment.h"

extern GeneralData* genData;

void BatchData::pre_filter_alignment_maps() {

    int num_removed = 0;
    int num_smplfy_reads = 0;
    bool removed_best = false;

    if (conf->verbose)
        fprintf(stdout, "\nPre-filtering alignment list ...\n");

    unordered_map<string, vector<Alignment> >::iterator r_idx;
    vector<Alignment>::iterator v_idx;
    set<vector<Alignment>::iterator> to_erase;

    for (r_idx = read_map_left.begin(); r_idx != read_map_left.end(); r_idx++) {
        if (r_idx->second.size() < 2) 
            continue;

        to_erase = filter_alignments(r_idx->second);
        num_removed += to_erase.size();

        removed_best = false;
        for (set<vector<Alignment>::iterator>::reverse_iterator e_idx = to_erase.rbegin(); e_idx != to_erase.rend(); e_idx++) {
            removed_best = (removed_best || (*e_idx)->is_best);
            r_idx->second.erase(*e_idx);
        }

        if (removed_best) {
            r_idx->second.begin()->is_best = true;
        }
        if (to_erase.size() > 0)
            num_smplfy_reads++;
    }

    to_erase.clear();
    for (r_idx = read_map_right.begin(); r_idx != read_map_right.end(); r_idx++) {
        if (r_idx->second.size() < 2) 
            continue;

        to_erase = filter_alignments(r_idx->second);
        num_removed += to_erase.size();

        removed_best = false;
        for (set<vector<Alignment>::iterator>::reverse_iterator e_idx = to_erase.rbegin(); e_idx != to_erase.rend(); e_idx++) {
            removed_best = (removed_best || (*e_idx)->is_best);
            r_idx->second.erase(*e_idx);
        }

        if (removed_best) {
            r_idx->second.begin()->is_best = true;
        }
        if (to_erase.size() > 0)
            num_smplfy_reads++;
    }
    if (conf->verbose)
        fprintf(stdout, "... removed %i alignments in %i reads.\n", num_removed, num_smplfy_reads);
}


void BatchData::parse_file() {

    char line[1000] ;
    char cp_line[1000];
    char* ret;

    unsigned int counter = 0;
    unsigned char pair_info = 0;
    
    Alignment curr_alignment;
    string id;
    string last_id = string("");

    FILE* infile = open_bam_pipe_in(conf->infile);
    bool read_header = false;
    bool unmapped = false;

    if (conf->verbose) {
        fprintf(stdout, "\nReading input file from: %s\n", conf->infile.c_str());
    }

    while (true) {

        ret = fgets(line, sizeof(line), infile);
        strcpy(cp_line, line);
        counter++;

        if (!ret) {
            if (counter == 1) {
                fprintf(stderr, "ERROR: Could not read SAM file %s\n", conf->infile.c_str());
                exit(1);
            }
            else {
                break;
            }
        }

        if (conf->verbose && counter % 1000000 == 0) {
            fprintf(stdout, "\n\t%i...", counter);
        }

        char* sl = strtok(line, "\t");
        
        string chr_name;

        if (line[0] == '@') {
            parse_header(sl);
            read_header = true;
            continue ;
        }

        if (! read_header) {
            fprintf(stderr, "\nERROR: Input file does not contain any header information!\n") ;
            exit(1);
        }
        
        unmapped = false;
        curr_alignment.clear();
        id = curr_alignment.fill(sl, pair_info, unmapped);

        if (unmapped)
            continue ;

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
            unordered_map<string, vector<Alignment> >::iterator id_idx = this->read_map_left.find(id);
            if (id_idx != this->read_map_left.end()) {
                id_idx->second.push_back(curr_alignment);
            } else {
                // create new element in read map
                vector<Alignment> tmp_align;
                tmp_align.push_back(curr_alignment);
                this->read_map_left.insert( pair<string, vector<Alignment> >(id, tmp_align) );
            }
        } else {
            unordered_map<string, vector<Alignment> >::iterator id_idx = this->read_map_right.find(id);
            if (id_idx != this->read_map_right.end()) {
                id_idx->second.push_back(curr_alignment);
            } else {
                // create new element in read map
                vector<Alignment> tmp_align;
                tmp_align.push_back(curr_alignment);
                this->read_map_right.insert( pair<string, vector<Alignment> >(id, tmp_align) );
            }
        }
    }

    if (conf->verbose) 
        fprintf(stdout, "\nsuccessfully parsed %i lines\n", counter);

    fclose(infile);
}

void BatchData::get_active_read_set() {

    //double insert1 = 0.0;
    //double insert2 = 0.0;

    unordered_map<string, vector<Alignment> >::iterator r_idx;
    unordered_map<string, vector<Alignment> >::iterator l_idx;

    // if pair processing, identify all compatible pairs
    // pairs are compatible, if they show same chr, opposite strands and
    // have an inner distance within the insert size range
    // store all left reads that do not have a right mapping in left_singles
    for (l_idx = this->read_map_left.begin(); l_idx != this->read_map_left.end(); l_idx++) {
        // check, if left read has corresponding right read
        r_idx = this->read_map_right.find(l_idx->first);
        if (conf->use_pair_info && r_idx != this->read_map_right.end()) {
            if (l_idx->second.size() < 2 && r_idx->second.size() < 2)
                continue;
            bool pair_vecs_empty = true;
            bool best_pair_left = false;
            bool best_pair_right = false;
            for (vector<Alignment>::iterator lv_idx = l_idx->second.begin(); lv_idx != l_idx->second.end(); lv_idx++) {
                for (vector<Alignment>::iterator rv_idx = r_idx->second.begin(); rv_idx != r_idx->second.end(); rv_idx++) {
                    if ((rv_idx->chr == lv_idx->chr) && (rv_idx->reversed != lv_idx->reversed)) {
                        //insert1 = abs((double) lv_idx->get_end() - (double) rv_idx->start);
                        //insert2 = abs((double) rv_idx->get_end() - (double) lv_idx->start);
                        //if (insert1 <= (conf->insert_size * (1.0 + conf->insert_dev)) || insert2 <= (conf->insert_size * (1.0 + conf->insert_dev))) {
                            if (pair_vecs_empty) {
                                vector<vector<Alignment>::iterator> tmp_vec;
                                tmp_vec.push_back(lv_idx);
                                this->active_left_pair.insert(pair<string, vector<vector<Alignment>::iterator> >(l_idx->first, tmp_vec));
                                tmp_vec.clear();
                                tmp_vec.push_back(rv_idx);
                                this->active_right_pair.insert(pair<string, vector<vector<Alignment>::iterator> >(r_idx->first, tmp_vec));
                                pair_vecs_empty = false;
                            } else {
                                this->active_left_pair[l_idx->first].push_back(lv_idx);
                                this->active_right_pair[r_idx->first].push_back(rv_idx);
                            }
                            if (! best_pair_left)
                                best_pair_left = lv_idx->is_best;
                            if (! best_pair_right)
                                best_pair_right = rv_idx->is_best;
                        //}
                    }
                }
            }
            if (pair_vecs_empty) {
                if (l_idx->second.size() > 1)
                    this->active_left_single.push_back(l_idx);
                if (r_idx->second.size() > 1)
                    this->active_right_single.push_back(r_idx);
            } else {
                // check, if current best alignment is part of active left pairs
                if (! best_pair_left) {
                    bool broken = false;
                    for (vector<Alignment>::iterator lv_idx = l_idx->second.begin(); lv_idx != l_idx->second.end(); lv_idx++) {
                        if (lv_idx->is_best) {
                            lv_idx->is_best = false;
                            lv_idx->update_coverage_map(0);
                            broken = true;
                            break;
                        }
                    }
                    assert(broken);
                    this->active_left_pair[l_idx->first].front()->is_best = true;
                    this->active_left_pair[l_idx->first].front()->update_coverage_map(1);
                }
                // check, if current best alignment is part of active right pairs
                if (! best_pair_right) {
                    bool broken = false;
                    for (vector<Alignment>::iterator rv_idx = r_idx->second.begin(); rv_idx != r_idx->second.end(); rv_idx++) {
                        if (rv_idx->is_best) {
                            rv_idx->is_best = false;
                            rv_idx->update_coverage_map(0);
                            broken = true;
                            break;
                        }
                    }
                    assert(broken);
                    this->active_right_pair[r_idx->first].front()->is_best = true;
                    this->active_right_pair[r_idx->first].front()->update_coverage_map(1);
                }
            }
        } else {
           if (l_idx->second.size() > 1)
               this->active_left_single.push_back(l_idx); 
        }
    }
    // check all (remaining) right reads
    for (r_idx = this->read_map_right.begin(); r_idx != this->read_map_right.end(); r_idx++) {
        if (!conf->use_pair_info) { 
            if (r_idx->second.size() > 1)
                this->active_right_single.push_back(r_idx);
        } else {
            l_idx = this->read_map_left.find(r_idx->first);
            if (l_idx == this->read_map_left.end() && r_idx->second.size() > 1)  
                this->active_right_single.push_back(r_idx);
        }
    }
}

double BatchData::get_total_min_loss() {

    double sum_min_loss = 0.0; 

    for (unordered_map<string, vector<Alignment> >::iterator r_idx = this->read_map_left.begin(); r_idx != this->read_map_left.end(); r_idx++) {
        for (vector<Alignment>::iterator v_idx = r_idx->second.begin(); v_idx != r_idx->second.end(); v_idx++) {
            if (v_idx->is_best) {
                vector<vector<unsigned long> > cov_keep;
                vector<vector<unsigned long> > cov_change;
                vector<pair<vector<Alignment>::iterator, bool> > aligns;
                aligns.push_back( make_pair(v_idx, true) );
                compute_coverage_loss(aligns, cov_keep, cov_change);
                sum_min_loss += get_variance(cov_keep);
            }
        }
    }
    for (unordered_map<string, vector<Alignment> >::iterator r_idx = this->read_map_right.begin(); r_idx != this->read_map_right.end(); r_idx++) {
        for (vector<Alignment>::iterator v_idx = r_idx->second.begin(); v_idx != r_idx->second.end(); v_idx++) {
            if (v_idx->is_best) {
                vector<vector<unsigned long>>  cov_keep;
                vector<vector<unsigned long> > cov_change;
                vector<pair<vector<Alignment>::iterator, bool> > aligns;
                aligns.push_back( make_pair(v_idx, true) );
                compute_coverage_loss(aligns, cov_keep, cov_change);
                sum_min_loss += get_variance(cov_keep);
            }
        }
    }
 
    if(conf->verbose)
        fprintf(stderr, "\n");

    return sum_min_loss;
}

unsigned int BatchData::smooth_coverage_map_single(unsigned int &num_ambiguous) {
    unsigned int num_changed = 0;
    num_changed += smooth_coverage_map_single_wrapper(this->active_left_single, num_ambiguous);
    num_changed += smooth_coverage_map_single_wrapper(this->active_right_single, num_ambiguous);
    return num_changed;
}

unsigned int BatchData::smooth_coverage_map_single_wrapper(list<unordered_map <string, vector<Alignment> >::iterator > &active_reads, unsigned int &num_ambiguous) {

        unsigned int num_changed = 0;
        unsigned int num_best = 0;
        double loss = 0.0;

        list<unordered_map<string, vector<Alignment> >::iterator>::iterator r_idx;
        for (r_idx = active_reads.begin(); r_idx != active_reads.end(); r_idx++) {
            if ((*r_idx)->second.size() == 1) {
                num_best++;
                continue;
            }

            num_ambiguous++;
            vector<Alignment>::iterator curr_best;
            vector<Alignment>::iterator v_idx;
            for (v_idx = (*r_idx)->second.begin(); v_idx != (*r_idx)->second.end(); v_idx++) {
                if (v_idx->is_best) {
                    curr_best = v_idx;
                    num_best++;
                }
            }

            bool changed = false;
            for (v_idx = (*r_idx)->second.begin(); v_idx != (*r_idx)->second.end(); v_idx++) {
                if (v_idx == curr_best)
                    continue;

                // check if v_idx < curr_best
                if (compare_single(v_idx, curr_best, loss)) {
                    changed = true;
                    curr_best->is_best = false;
                    curr_best->update_coverage_map(0);
                    curr_best = v_idx;          
                    curr_best->is_best = true;
                    curr_best->update_coverage_map(1);
                }
            }
            if (changed)
                num_changed++;
        }

       if ((size_t) num_best != active_reads.size()) {
            fprintf(stderr, "best: %i  map size: %i", num_best, (int) active_reads.size());
            assert((size_t) num_best == active_reads.size());
        }

        return num_changed;
}

unsigned int BatchData::smooth_coverage_map_paired(unsigned int &num_ambiguous) {

        unsigned int num_changed = 0;
        double loss = 0.0;

        unordered_map<string, vector<vector<Alignment>::iterator> >::iterator l_idx;
        unordered_map<string, vector<vector<Alignment>::iterator> >::iterator r_idx;

        vector<vector<Alignment>::iterator>::iterator lv_idx;
        vector<vector<Alignment>::iterator>::iterator rv_idx;

        vector<Alignment>::iterator best_left_idx;
        vector<Alignment>::iterator best_right_idx;
        bool found_best = false;
        bool changed = false;

        for (l_idx = this->active_left_pair.begin(), r_idx = this->active_right_pair.begin(); l_idx != this->active_left_pair.end() && r_idx != this->active_right_pair.end(); l_idx++, r_idx++) {
            // find a best pairing
            found_best = false;
            changed = false;
            // there is only one pair
            if (l_idx->second.size()==1) {
                continue;
            }
            num_ambiguous++;
            for(lv_idx = l_idx->second.begin(), rv_idx = r_idx->second.begin(); lv_idx < l_idx->second.end() && rv_idx  < r_idx->second.end(); lv_idx++, rv_idx++) {
                if ((*rv_idx)->is_best && (*lv_idx)->is_best) {
                    best_left_idx = (*lv_idx);
                    best_right_idx = (*rv_idx);
                    found_best = true;
                    break;
                }
            }
            if (!found_best) {
                for(lv_idx = l_idx->second.begin(), rv_idx = r_idx->second.begin(); lv_idx < l_idx->second.end() && rv_idx  < r_idx->second.end(); lv_idx++, rv_idx++) {
                    (*lv_idx)->is_best = false;
                    (*rv_idx)->is_best = false;
                }
                best_left_idx = l_idx->second.front();
                best_right_idx = r_idx->second.front();

                best_left_idx->is_best = true;
                best_right_idx->is_best = true;
            }
            for(lv_idx = l_idx->second.begin(), rv_idx = r_idx->second.begin(); lv_idx < l_idx->second.end() && rv_idx  < r_idx->second.end(); lv_idx++, rv_idx++) {
                if ((*lv_idx)  == best_left_idx && (*rv_idx) == best_right_idx)
                    continue;
                if (compare_pair(*lv_idx, *rv_idx, best_left_idx, best_right_idx, loss)) {
                    changed = true;
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
            }

            if (changed)
                num_changed++;
        
        }
      /*  map <int, vector<unsigned short> >::iterator it = genData->coverage_map.begin();
        for (it; it != genData->coverage_map.end(); it++) {
            fprintf(stdout, "cov vec:\n");
            for (vector<unsigned short>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                fprintf(stdout, "%i ", (*it2));
            }
            fprintf(stdout, "\n");
        }*/

        /*if ((size_t) num_best != read_map.size()) {
            fprintf(stderr, "best: %i  map size: %i", num_best, (int) read_map.size());
            assert((size_t) num_best == read_map.size());
        }*/
        //fprintf(stdout, "changed %i in paired", num_changed);
        return num_changed;
}

