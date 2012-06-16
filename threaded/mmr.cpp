#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <assert.h>
#include <pthread.h>

#include <unordered_map>
#include <map>
#include <vector>
#include <list>
#include <set>
#include <queue>

#include "config.h"
#include "alignment.h"
#include "data.h"

using namespace std;

// global variables
GeneralData* genData = new GeneralData::GeneralData();
Config* conf;
bool done;

pthread_mutex_t mutex_coverage;
pthread_mutex_t mutex_fifo;
pthread_mutex_t mutex_done;
pthread_mutex_t mutex_best_left;
pthread_mutex_t mutex_best_right;

pthread_cond_t fifo_is_full;
pthread_attr_t attr;

queue<ThreadData_Online*> fifo;

typedef struct {

} batch_data;

typedef struct {
    
} online_data;

FILE* open_bam_pipe_in(std::string & in_fname)
{
	string command = string("/fml/ag-raetsch/share/software/samtools-0.1.7a/samtools view -h 2> /dev/null ") + in_fname  + " && echo samtools subprocess for reading terminated successfully";
	FILE* IN_FP=NULL ;
	fflush(stdout) ;
	IN_FP = popen(command.c_str(), "r") ;
	return IN_FP ;
}

FILE* open_bam_pipe_out(std::string & out_fname)
{
	string command = string("/fml/ag-raetsch/share/software/samtools-0.1.7a/samtools view -bS -o ") + out_fname +  " - 2> /dev/null" + " && echo samtools for writing subprocess terminated successfully";
	FILE* OUT_FP=NULL ;
	fflush(stdin) ;
	OUT_FP = popen(command.c_str(), "w") ;
	return OUT_FP ;
}

void parse_cigar(string cigar, vector<char> &operations, vector<int> &sizes) {
    
    size_t op_pos = cigar.find_first_of("HSMIDNXP");
    while (op_pos != string::npos) {
        sizes.push_back(atoi(cigar.substr(0, op_pos).c_str()));
        operations.push_back(cigar.at(op_pos));
        cigar = cigar.substr(op_pos + 1, cigar.size());
        op_pos = cigar.find_first_of("HSMIDNXP");
    }
}

unsigned int get_alignment_end(string cigar, unsigned int start) {
    
    unsigned int end = start;
    size_t op_pos = cigar.find_first_of("HSMIDNXP");
    while (op_pos != string::npos) {
        if ((cigar.at(op_pos) == 'N') || (cigar.at(op_pos) == 'D') || (cigar.at(op_pos) == 'M')) {
            end += (atoi(cigar.substr(0, op_pos).c_str()));
        }
        cigar = cigar.substr(op_pos + 1, cigar.size());
        op_pos = cigar.find_first_of("HSMIDNXP");
    }
    return end;
}

void update_coverage_map(const vector<alignment>::iterator a_it, vector<unsigned short> &coverage_map, bool positive) {

    pthread_mutex_lock(&mutex_coverage);
    vector<unsigned short>::iterator idx = coverage_map.begin() + a_it->start;

    vector<char> operations;
    vector<int> sizes;
    
    parse_cigar(a_it->cigar, operations, sizes);

    for (size_t i = 0; i < sizes.size(); i++) {
        switch (operations.at(i)) {
            case 'M': case 'D': for (int j = 0; j < sizes.at(i); j++) {if (idx < coverage_map.end()) {*idx += (*idx > 0 || positive)?(2*positive - 1):0; idx++;}}; break;
            case 'N': idx += sizes.at(i);
        }
        if (idx >= coverage_map.end())
            break;
    }
    pthread_mutex_unlock(&mutex_coverage);
}

vector<unsigned short> get_coverage(vector<unsigned short> &coverage_map, string cigar, unsigned long start, unsigned int window_size, vector<unsigned short> &intron_cov) {

    vector<char> operations;
    vector<int> sizes;

    vector<unsigned short> curr_cov;
    vector<unsigned short>::iterator cov_idx = coverage_map.begin();
    vector<unsigned short>::iterator internal_cov_idx;

    unsigned int offset = (window_size>=start)?start:window_size;
    size_t step_size = 1;

    if (window_size < start) 
        cov_idx += (start - window_size);
        
    parse_cigar(cigar, operations, sizes);

    for (size_t i = 0; i < offset; i++) {
        if (cov_idx < coverage_map.end()) {
            curr_cov.push_back(*cov_idx);
            cov_idx++;
        }
        else {
            return curr_cov;
        }
    }

    for (size_t i = 0; i < sizes.size(); i++) {
        switch (operations.at(i)) {
            case 'M': case 'D': for (int j = 0; j < sizes.at(i); j++) {if (cov_idx < coverage_map.end()) {curr_cov.push_back(*cov_idx); cov_idx++;}}; break;
            case 'N': { step_size = max(sizes.at(i) / 50, 1);
                        internal_cov_idx = cov_idx + conf->intron_offset;
                        for (int j = conf->intron_offset; j < sizes.at(i) - conf->intron_offset; j+=step_size) {
                            if (internal_cov_idx < coverage_map.end()) {
                                intron_cov.push_back(*internal_cov_idx); 
                                internal_cov_idx++;
                            }
                        }; 
                        cov_idx += sizes.at(i);
                        break;
                      }
        }
    }

    for (size_t i = 0; i < window_size; i++) {
        if (cov_idx < coverage_map.end()) {
            curr_cov.push_back(*cov_idx);
            cov_idx++;
        }
        else {
            break;
        }
    }
    return curr_cov;
}

double intron_penalty(vector<unsigned short> &intron_coverage) {
    double sum = 0.0;

    if (intron_coverage.size() > 0) {
        for (vector<unsigned short>::iterator it = intron_coverage.begin(); it != intron_coverage.end(); it++) {
            sum += (double) *it;
        }
        sum /= (double) intron_coverage.size();
    }

    return sum;
}

double get_variance(vector<unsigned short> &coverage, vector<unsigned short> &intron_coverage) {

    if (coverage.size() < 2) {
        return 0.0;
    } else {
        double sum = 0.0;
        double int_pnlty = 0.0;
        for ( vector<unsigned short>::iterator it = coverage.begin(); it != coverage.end(); it++) {
            sum += (double) *it;
        }
        double mean = sum / (double) coverage.size();
        sum = 0.0;
        for ( vector<unsigned short>::iterator it = coverage.begin(); it != coverage.end(); it++) {
            sum += pow((double) *it - mean, 2.0);
        }

        if (intron_coverage.size() > 0) {
            int_pnlty = intron_penalty(intron_coverage);
        }

        return (sum / ((double) coverage.size() - 1.0)) + int_pnlty;  
    }
}

vector<unsigned short> alter_coverage(vector<unsigned short> &source, unsigned int window_left, unsigned int window_right, bool is_positive) {

    vector<unsigned short> target;

    for (size_t i = 0; i < source.size(); i++) {
        if (i < window_left || i >= source.size() - window_right) {
            target.push_back(source.at(i));
        } else {
            target.push_back((source.at(i) > 0 || is_positive)?source.at(i) + (2*is_positive) - 1:0);
        }
    }
    return target;
}

bool compare_pair(vector<alignment>::iterator candidate_left, vector<alignment>::iterator candidate_right, vector<alignment>::iterator best_left, vector<alignment>::iterator best_right, map <int, vector<unsigned short> > &coverage_map, Config* conf) {

    vector<unsigned short> intron_cov_candidate_left;
    vector<unsigned short> intron_cov_candidate_right;
    vector<unsigned short> intron_cov_best_left;
    vector<unsigned short> intron_cov_best_right;

    vector<unsigned short> exon_cov_candidate_left_without = get_coverage(coverage_map[candidate_left->chr], candidate_left->cigar, candidate_left->start, conf->window_size, intron_cov_candidate_left);
    vector<unsigned short> exon_cov_candidate_left_with = alter_coverage(exon_cov_candidate_left_without, min((unsigned int) candidate_left->start, conf->window_size), conf->window_size, true);
    vector<unsigned short> exon_cov_candidate_right_without = get_coverage(coverage_map[candidate_right->chr], candidate_right->cigar, candidate_right->start, conf->window_size, intron_cov_candidate_right);
    vector<unsigned short> exon_cov_candidate_right_with = alter_coverage(exon_cov_candidate_right_without, min((unsigned int) candidate_right->start, conf->window_size), conf->window_size, true);

    vector<unsigned short> exon_cov_best_left_with = get_coverage(coverage_map[best_left->chr], best_left->cigar, best_left->start, conf->window_size, intron_cov_best_left);
    vector<unsigned short> exon_cov_best_left_without = alter_coverage(exon_cov_best_left_with, min((unsigned int) best_left->start, conf->window_size), conf->window_size, false);
    vector<unsigned short> exon_cov_best_right_with = get_coverage(coverage_map[best_right->chr], best_right->cigar, best_right->start, conf->window_size, intron_cov_best_right);
    vector<unsigned short> exon_cov_best_right_without = alter_coverage(exon_cov_best_right_with, min((unsigned int) best_right->start, conf->window_size), conf->window_size, false);

    vector<unsigned short> empty_cov;
    double var_cand_left_without = get_variance(exon_cov_candidate_left_without, empty_cov);
    double var_cand_left_with = get_variance(exon_cov_candidate_left_with, intron_cov_candidate_left);
    double var_cand_right_without = get_variance(exon_cov_candidate_right_without, empty_cov);
    double var_cand_right_with = get_variance(exon_cov_candidate_right_with, intron_cov_candidate_right);
    
    double var_best_left_without = get_variance(exon_cov_best_left_without, empty_cov);
    double var_best_left_with = get_variance(exon_cov_best_left_with, intron_cov_best_left);
    double var_best_right_without = get_variance(exon_cov_best_right_without, empty_cov);
    double var_best_right_with = get_variance(exon_cov_best_right_with, intron_cov_best_right);

    return (var_cand_left_with + var_cand_right_with + var_best_left_without + var_best_right_without) < (var_cand_left_without + var_cand_right_without + var_best_left_with + var_best_right_with);
}

bool compare_single(vector<alignment>::iterator candidate, vector<alignment>::iterator best, map <int, vector<unsigned short> > &coverage_map, Config* conf) {

    vector<unsigned short> intron_cov_candidate;
    vector<unsigned short> intron_cov_best;
    vector<unsigned short> exon_cov_candidate_without = get_coverage(coverage_map[candidate->chr], candidate->cigar, candidate->start, conf->window_size, intron_cov_candidate);
    vector<unsigned short> exon_cov_best_with = get_coverage(coverage_map[best->chr], best->cigar, best->start, conf->window_size, intron_cov_best);
    vector<unsigned short> exon_cov_candidate_with = alter_coverage(exon_cov_candidate_without, min((unsigned int) candidate->start, conf->window_size), conf->window_size, true);
    vector<unsigned short> exon_cov_best_without = alter_coverage(exon_cov_best_with, min((unsigned int) best->start, conf->window_size), conf->window_size, false);

    vector<unsigned short> empty_cov;
    double var_cand_without = get_variance(exon_cov_candidate_without, empty_cov);
    double var_best_with = get_variance(exon_cov_best_with, intron_cov_best);
    double var_cand_with = get_variance(exon_cov_candidate_with, intron_cov_candidate);
    double var_best_without = get_variance(exon_cov_best_without, empty_cov);

    return (var_cand_with + var_best_without) < (var_cand_without + var_best_with);
}

void parse_header(char* sl) {

    int idx = 0;
    string chr_name;
    while(sl != NULL) {
        string tmp_sl = string(sl);
        
        if (idx == 1) {
            tmp_sl = tmp_sl.substr(3, tmp_sl.size());
            chr_name = tmp_sl ;
            if (genData->chr_num.find(tmp_sl) == genData->chr_num.end()) {
                genData->chr_num.insert( pair<string, int>(tmp_sl,  genData->chr_num.size() + 1) );
            }
            else {
                fprintf(stderr, "WARNING: Doubled contig names in header!\n Ignoring %s\n\n", tmp_sl.c_str());
            }
        }
        else if (idx == 2) {
            tmp_sl = tmp_sl.substr(3, tmp_sl.size());
            genData->chr_size.push_back(atoi(tmp_sl.c_str()));    
            vector<unsigned short> tmp_cov(genData->chr_size.back(), 0);
            if (conf->verbose) 
                fprintf(stdout, "\t...reserving memory for contig %s of size %i\n", chr_name.c_str(), genData->chr_size.back());
            genData->coverage_map.insert( pair<int, vector<unsigned short> >(genData->chr_num[chr_name], tmp_cov) ); 
        }   
        idx ++;
        sl = strtok(NULL, "\t");    
    }

    if (genData->chr_size.size() != genData->chr_num.size() || genData->chr_num.size() != genData->coverage_map.size()) {
        fprintf(stderr, "\nERROR: Header information incomplete!");
        exit(-1);
    }
}

string get_alignment(char* sl, alignment &curr_alignment, map<string, int> &chr_num, unsigned char &pair) {

    int chr_n = 0;
    int idx = 0;
    unsigned long start = 0;
    unsigned char edit_ops = 0;
    unsigned char quality = 0;
    string id, chr, cigar;
    bool reversed = false;

    while (sl != NULL) {
        if (idx == 0) { 
            id = sl;
        } else if (idx == 1) {
            pair = (atoi(sl) & 128);
            reversed = ((atoi(sl) & 16) == 16);
        } else if (idx == 2) {
            if (chr_num.find(sl) == chr_num.end()) {
                fprintf(stderr, "ERROR: Contig name not in header!\n Contig: %s\n\n", sl);
                exit(1);
            } else {
                chr_n = chr_num[sl];
            }
        } else if (idx == 3) {
            start = strtoul(sl, NULL, 0) - 1;
        } else if (idx == 4) {
            quality = (unsigned  char) atoi(sl);
        } else if (idx == 5) {
            cigar = sl;
        } else if (!conf->use_variants && strlen(sl) > 3 && sl[0] == 'N' && sl[1] == 'M' && sl[2] == ':') {
            edit_ops = atoi(sl+5); 
            break;
        } else if (conf->use_variants && strlen(sl) > 3 && sl[0] == 'X' && (sl[1] == 'M' || sl[1] == 'G') && sl[2] == ':') {
            edit_ops += atoi(sl+5);
            break;
        }
        sl = strtok(NULL, "\t");
        idx ++;
    }
    if (idx < 6) {
        return(string(""));
    }
    curr_alignment = alignment(chr_n, start, cigar, false, edit_ops, quality, reversed);
    return id;
}

string update_line_flag(char* line, bool is_best) {

    char cp_line[1000];
    strcpy(cp_line, line);

    char* sl = strtok(line, "\t");

    int idx = 0;
    string return_line;
    while (sl != NULL) {
        if (idx == 1) {
            if ( (!is_best && (atoi(sl) & 256) == 256) || (is_best && (atoi(sl) & 256) == 0)) { 
                return_line = string(cp_line);
                return return_line.substr(0, return_line.size() - 1);
            }
            else {
                char new_flag[10]; 
                if (is_best) {
                    sprintf(new_flag, "%i", atoi(sl) - 256);        
                }
                else {
                    sprintf(new_flag, "%i", atoi(sl) + 256);        
                }
                return_line += (string("\t") + new_flag);
            }
        }
        else {
            return_line += (string("\t") + sl);
        }
        sl = strtok(NULL, "\t");
        idx++;
    }
    return return_line.substr(1, return_line.size());
}

set<vector<alignment>::iterator> filter_alignments(vector<alignment> &aligns) {

    set<vector<alignment>::iterator> to_erase;
    vector<alignment>::iterator v_idx;
    //sort(aligns.begin(), aligns.end(), alignment_comparator);
    if (aligns.size() == 0)
        return to_erase;
    unsigned int min_ops = aligns.begin()->edit_ops;
    for (v_idx = aligns.begin() + 1; v_idx < aligns.end(); v_idx++)
        min_ops = v_idx->edit_ops < min_ops?v_idx->edit_ops:min_ops;

    for (v_idx = aligns.begin(); v_idx < aligns.end(); v_idx++) {
        if (v_idx->edit_ops - min_ops > conf->filter_distance) {
            to_erase.insert(v_idx);
        }
    }
    return to_erase;
}

void pre_filter_alignment_map(unordered_map<string, vector<alignment> > &read_map) {

    int num_removed = 0;
    int num_smplfy_reads = 0;
    bool removed_best = false;

    if (conf->verbose)
        fprintf(stdout, "\nPre-filtering alignment list ...\n");

    unordered_map<string, vector<alignment> >::iterator r_idx;
    vector<alignment>::iterator v_idx;
    set<vector<alignment>::iterator> to_erase;

    for (r_idx = read_map.begin(); r_idx != read_map.end(); r_idx++) {
        if (r_idx->second.size() < 2) 
            continue;

        to_erase = filter_alignments(r_idx->second);
        num_removed += to_erase.size();

        removed_best = false;
        for (set<vector<alignment>::iterator>::reverse_iterator e_idx = to_erase.rbegin(); e_idx != to_erase.rend(); e_idx++) {
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

void write_output_direct() {
//unordered_map <string, size_t, hash<string> > &best_left, unordered_map <string, size_t, hash<string> > &best_right, map<string, int> &chr_num) {

    FILE* outfile = open_bam_pipe_out(conf->outfile);
    FILE* infile = open_bam_pipe_in(conf->infile);

    char* ret;
    char line[1000];
    char cp_line[1000];

    unsigned char pair_info = 0;

    size_t left_counter = 0;
    size_t right_counter = 0;

    int output_counter = 0;

    string id;
    string last_left_id = string("");
    string last_right_id = string("");
    alignment curr_alignment;

    if (conf->outfile.size() == 0) {
        fprintf(stderr, "\nERROR: No outfile defined!\n\n");
        exit(-1);
    }

    if (conf->verbose) { 
        fprintf(stdout, "Writing output to %s ...\n", conf->outfile.c_str());
    }
    
    while (true) {

        ret = fgets(line, sizeof(line), infile);
        if (!ret)
            break;
        strcpy(cp_line, line);

        if (line[0] == '@') {
            fprintf(outfile, "%s", line);
            continue;
        }

        char* sl = strtok(line, "\t");
        id = get_alignment(sl, curr_alignment, genData->chr_num, pair_info);
        if (id.size() == 0) {
            continue ;
        }
    
        if (pair_info == 0) {
            if (id.compare(last_left_id))
                left_counter = 1;
            else
                left_counter++;
            last_left_id = id;

            if (genData->best_left.find(id) == genData->best_left.end())
                continue;
            else {
                if (!conf->print_best_only || (genData->best_left[id] == (left_counter - 1))) { 
                    string print_line = update_line_flag(cp_line, (genData->best_left[id] == (left_counter - 1)));
                    fprintf(outfile, "%s\n", print_line.c_str());
                    output_counter++;
                }
            }
        } else {
            if (id.compare(last_right_id))
                right_counter = 1;
            else
                right_counter++;
            last_right_id = id;

            if (genData->best_right.find(id) == genData->best_right.end())
                continue;
            else {
                if (!conf->print_best_only || (genData->best_right[id] == (right_counter - 1))) { 
                    string print_line = update_line_flag(cp_line, (genData->best_right[id] == (right_counter - 1)));
                    fprintf(outfile, "%s\n", print_line.c_str());
                    output_counter++;
                }
            }
        }

    }
    fclose(outfile);
    fclose(infile);

    if (conf->verbose) 
        fprintf(stdout, "... done.\nPrinted %i lines.\n", output_counter);
}


void write_output(unordered_map <string, vector<alignment> > &read_map_left, unordered_map <string, vector<alignment> > &read_map_right, map<string, int> &chr_num) {

    FILE* outfile = open_bam_pipe_out(conf->outfile);
    FILE* infile = open_bam_pipe_in(conf->infile);

    char* ret;
    char line[1000];
    char cp_line[1000];

    unsigned char pair_info = 0;

    int output_counter = 0;

    string id;
    alignment curr_alignment;

    if (conf->outfile.size() == 0) {
        fprintf(stderr, "\nERROR: No outfile defined!\n\n");
        exit(-1);
    }

    if (conf->verbose) { 
        fprintf(stdout, "Writing output to %s ...\n", conf->outfile.c_str());
    }

    while (true) {

        ret = fgets(line, sizeof(line), infile);
        if (!ret)
            break;
        strcpy(cp_line, line);

        if (line[0] == '@') {
            fprintf(outfile, "%s", line);
            continue;
        }

        char* sl = strtok(line, "\t");
        id = get_alignment(sl, curr_alignment, chr_num, pair_info);
        if (id.size() == 0) {
            continue ;
        }
    
        vector<alignment>::iterator it;
        if (pair_info == 0) {
            it = find(read_map_left[id].begin(), read_map_left[id].end(), curr_alignment) ;
            if (it == read_map_left[id].end())
                continue;
        } else {
            it = find(read_map_right[id].begin(), read_map_right[id].end(), curr_alignment) ;
            if (it == read_map_right[id].end())
                continue;
        }

        if (!conf->print_best_only || it->is_best) { 
            string print_line = update_line_flag(cp_line, it->is_best);
            fprintf(outfile, "%s\n", print_line.c_str());
            output_counter++;
        }
    }
    fclose(outfile);
    fclose(infile);

    if (conf->verbose) 
        fprintf(stdout, "... done.\nPrinted %i lines.\n", output_counter);
}

double get_total_min_loss(unordered_map<string, vector<alignment> > &read_map, map <int, vector<unsigned short> > &coverage_map) {

    double sum_min_loss = 0.0; 
    vector<unsigned short> tmp_cov;
    vector<unsigned short> intron_cov;

    for (unordered_map<string, vector<alignment> >::iterator r_idx = read_map.begin(); r_idx != read_map.end(); r_idx++) {
        for (vector<alignment>::iterator v_idx = r_idx->second.begin(); v_idx != r_idx->second.end(); v_idx++) {
            if (v_idx->is_best) {
                intron_cov.clear();
                tmp_cov = get_coverage(coverage_map[v_idx->chr], v_idx->cigar, v_idx->start, conf->window_size, intron_cov);
                sum_min_loss += get_variance(tmp_cov, intron_cov); 
            }
        }
    }
    if(conf->verbose)
        fprintf(stderr, "\n");

    return sum_min_loss;

}

unsigned int smooth_coverage_map_single(list<unordered_map <string, vector<alignment> >::iterator > &active_reads, map<int, vector<unsigned short> > &coverage_map, unsigned int &num_ambiguous) {

        unsigned int num_changed = 0;
        unsigned int num_best = 0;

        list<unordered_map<string, vector<alignment> >::iterator>::iterator r_idx;
        for (r_idx = active_reads.begin(); r_idx != active_reads.end(); r_idx++) {
            if ((*r_idx)->second.size() == 1) {
                num_best++;
                continue;
            }

            num_ambiguous++;
            vector<alignment>::iterator curr_best;
            vector<alignment>::iterator v_idx;
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
                if (compare_single(v_idx, curr_best, coverage_map, conf)) {
                    changed = true;
                    curr_best->is_best = false;
                    update_coverage_map(curr_best, coverage_map[curr_best->chr], 0);
                    curr_best = v_idx;          
                    curr_best->is_best = true;
                    update_coverage_map(curr_best, coverage_map[curr_best->chr], 1);
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

unsigned int smooth_coverage_map_paired(unordered_map <string, vector<vector<alignment>::iterator> > &active_left_pair, unordered_map <string, vector<vector<alignment>::iterator> > &active_right_pair, map<int, vector<unsigned short> > &coverage_map, unsigned int &num_ambiguous) {

        unsigned int num_changed = 0;
        //unsigned int num_best = 0;

        unordered_map<string, vector<vector<alignment>::iterator> >::iterator l_idx;
        unordered_map<string, vector<vector<alignment>::iterator> >::iterator r_idx;

        vector<vector<alignment>::iterator>::iterator lv_idx;
        vector<vector<alignment>::iterator>::iterator rv_idx;

        vector<alignment>::iterator best_left_idx;
        vector<alignment>::iterator best_right_idx;
        bool found_best = false;
        bool changed = false;

        for (l_idx = active_left_pair.begin(), r_idx = active_left_pair.begin(); l_idx != active_left_pair.end() && r_idx != active_right_pair.end(); l_idx++, r_idx++) {
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
                if (compare_pair(*lv_idx, *rv_idx, best_left_idx, best_right_idx, coverage_map, conf)) {
                    changed = true;
                    best_left_idx->is_best = false;
                    best_right_idx->is_best = false;
                    update_coverage_map(best_left_idx, coverage_map[best_left_idx->chr], 0);
                    update_coverage_map(best_right_idx, coverage_map[best_left_idx->chr], 0);
                    best_left_idx = *lv_idx;
                    best_right_idx = *rv_idx;
                    best_left_idx->is_best = true;
                    best_right_idx->is_best = true;
                    update_coverage_map(best_left_idx, coverage_map[best_left_idx->chr], 1);
                    update_coverage_map(best_right_idx, coverage_map[best_left_idx->chr], 1);
                }
            }

            if (changed)
                num_changed++;
        
        }

        /*if ((size_t) num_best != read_map.size()) {
            fprintf(stderr, "best: %i  map size: %i", num_best, (int) read_map.size());
            assert((size_t) num_best == read_map.size());
        }*/
        //fprintf(stdout, "changed %i in paired", num_changed);
        return num_changed;
}

void get_active_read_set(ThreadData_Batch* data) {

    double insert1 = 0.0;
    double insert2 = 0.0;

    unordered_map<string, vector<alignment> >::iterator r_idx;
    unordered_map<string, vector<alignment> >::iterator l_idx;

    // if pair processing, identify all compatible pairs
    // pairs are compatible, if they show same chr, opposite strands and
    // have an inner distance within the insert size range
    // store all left reads that do not have a right mapping in left_singles
    for (l_idx = data->read_map_left.begin(); l_idx != data->read_map_left.end(); l_idx++) {
        // check, if left read has corresponding right read
        r_idx = data->read_map_right.find(l_idx->first);
        if (conf->use_pair_info && r_idx != data->read_map_right.end()) {
            if (l_idx->second.size() < 2 && r_idx->second.size() < 2)
                continue;
            bool inserted = false;
            bool best_pair_left = false;
            bool best_pair_right = false;
            for (vector<alignment>::iterator lv_idx = l_idx->second.begin(); lv_idx != l_idx->second.end(); lv_idx++) {
                for (vector<alignment>::iterator rv_idx = r_idx->second.begin(); rv_idx != r_idx->second.end(); rv_idx++) {
                    if ((rv_idx->chr == lv_idx->chr) && (rv_idx->reversed != lv_idx->reversed)) {
                        insert1 = abs((double) get_alignment_end(lv_idx->cigar, lv_idx->start) - (double) rv_idx->start);
                        insert2 = abs((double) get_alignment_end(rv_idx->cigar, rv_idx->start) - (double) lv_idx->start);
                        if (insert1 <= (conf->insert_size * (1.0 + conf->insert_dev)) || insert2 <= (conf->insert_size * (1.0 + conf->insert_dev))) {
                            if (!inserted) {
                                vector<vector<alignment>::iterator> tmp_vec;
                                tmp_vec.push_back(lv_idx);
                                data->active_left_pair.insert(pair<string, vector<vector<alignment>::iterator> >(l_idx->first, tmp_vec));
                                tmp_vec.clear();
                                tmp_vec.push_back(rv_idx);
                                data->active_right_pair.insert(pair<string, vector<vector<alignment>::iterator> >(r_idx->first, tmp_vec));
                                inserted = true;
                            } else {
                                data->active_left_pair[l_idx->first].push_back(lv_idx);
                                data->active_right_pair[r_idx->first].push_back(rv_idx);
                            }
                            if (! best_pair_left)
                                best_pair_left = lv_idx->is_best;
                            if (! best_pair_right)
                                best_pair_right = rv_idx->is_best;
                        }
                    }
                }
            }
            if (!inserted) {
                if (l_idx->second.size() > 1)
                    data->active_left_single.push_back(l_idx);
                if (r_idx->second.size() > 1)
                    data->active_right_single.push_back(r_idx);
            } else {
                // check, if current best alignment is part of active left pairs
                if (! best_pair_left) {
                    for (vector<alignment>::iterator lv_idx = l_idx->second.begin(); lv_idx != l_idx->second.end(); lv_idx++) {
                        if (lv_idx->is_best) {
                            lv_idx->is_best = false;
                            update_coverage_map(lv_idx, genData->coverage_map[lv_idx->chr], 0);
                            break;
                        }
                    }
                    data->active_left_pair[l_idx->first].front()->is_best = true;
                    update_coverage_map(data->active_left_pair[l_idx->first].front(), genData->coverage_map[data->active_left_pair[l_idx->first].front()->chr], 1);
                }
                // check, if current best alignment is part of active right pairs
                if (! best_pair_right) {
                    for (vector<alignment>::iterator rv_idx = r_idx->second.begin(); rv_idx != r_idx->second.end(); rv_idx++) {
                        if (rv_idx->is_best) {
                            rv_idx->is_best = false;
                            update_coverage_map(rv_idx, genData->coverage_map[rv_idx->chr], 0);
                            break;
                        }
                    }
                    data->active_right_pair[r_idx->first].front()->is_best = true;
                    update_coverage_map(data->active_right_pair[r_idx->first].front(), genData->coverage_map[data->active_right_pair[r_idx->first].front()->chr], 1);
                }
            }
        } else {
           if (l_idx->second.size() > 1)
               data->active_left_single.push_back(l_idx); 
        }
    }
    // check all (remaining) right reads
    for (r_idx = data->read_map_right.begin(); r_idx != data->read_map_right.end(); r_idx++) {
        if (!conf->use_pair_info) { 
            if (r_idx->second.size() > 1)
                data->active_right_single.push_back(r_idx);
        } else {
            l_idx = data->read_map_left.find(r_idx->first);
            if (l_idx == data->read_map_left.end() && r_idx->second.size() > 1)  
                data->active_right_single.push_back(r_idx);
        }
    }
}

void get_active_reads(string read_id, set<vector<alignment>::iterator> &ignore_reads_left, set<vector<alignment>::iterator> &ignore_reads_right, vector<alignment> &left_reads, vector<alignment> &right_reads, vector<vector<alignment>::iterator> &active_left_reads, vector<vector<alignment>::iterator> &active_right_reads, unordered_map <string, size_t, hash<string> > &best_left, unordered_map <string, size_t, hash<string> > &best_right, bool &found_pairs) {

    double insert1 = 0.0;
    double insert2 = 0.0;

    // if pair processing, identify all compatible pairs
    // pairs are compatible, if they show same chr, opposite strands and
    // have an inner distance within the insert size range

    active_left_reads.clear();
    active_right_reads.clear();

    if (conf->use_pair_info && ! (right_reads.size() == 0 || left_reads.size() == 0)) {
        bool best_pair = false;
        for (vector<alignment>::iterator lv_idx = left_reads.begin(); lv_idx != left_reads.end(); lv_idx++) {
            if (conf->pre_filter && (ignore_reads_left.find(lv_idx) != ignore_reads_left.end()))
                continue;
            for (vector<alignment>::iterator rv_idx = right_reads.begin(); rv_idx != right_reads.end(); rv_idx++) {
                if (conf->pre_filter && (ignore_reads_right.find(rv_idx) != ignore_reads_right.end()))
                    continue;
                if ((rv_idx->chr == lv_idx->chr) && (rv_idx->reversed != lv_idx->reversed)) {
                    insert1 = abs((double) get_alignment_end(lv_idx->cigar, lv_idx->start) - (double) rv_idx->start);
                    insert2 = abs((double) get_alignment_end(rv_idx->cigar, rv_idx->start) - (double) lv_idx->start);
                    if (insert1 <= (conf->insert_size * (1.0 + conf->insert_dev)) || insert2 <= (conf->insert_size * (1.0 + conf->insert_dev))) {
                        if (! best_pair) 
                            best_pair = (lv_idx->is_best && rv_idx->is_best);
                        active_left_reads.push_back(lv_idx);
                        active_right_reads.push_back(rv_idx);
                        found_pairs = true;
                    }
                }
            }
        }

        // no active pair is best alignment
        if (! best_pair && found_pairs) {
            for (vector<alignment>::iterator lv_idx = left_reads.begin(); lv_idx != left_reads.end(); lv_idx++) {
                if (lv_idx->is_best) {
                    lv_idx->is_best = false;
                    update_coverage_map(lv_idx, genData->coverage_map[lv_idx->chr], 0);
                    break;
                }
            }
            active_left_reads.front()->is_best = true;
            update_coverage_map(active_left_reads.front(), genData->coverage_map[active_left_reads.front()->chr], 1);

            pthread_mutex_lock(&mutex_best_left);
            best_left[read_id] = (active_left_reads.front() - left_reads.begin());
            pthread_mutex_unlock(&mutex_best_left);

            for (vector<alignment>::iterator rv_idx = right_reads.begin(); rv_idx != right_reads.end(); rv_idx++) {
                if (rv_idx->is_best) {
                    rv_idx->is_best = false;
                    update_coverage_map(rv_idx, genData->coverage_map[rv_idx->chr], 0);
                    break;
                }
            }
            active_right_reads.front()->is_best = true;
            update_coverage_map(active_right_reads.front(), genData->coverage_map[active_right_reads.front()->chr], 1);

            pthread_mutex_lock(&mutex_best_right);
            assert(active_right_reads.front() - right_reads.begin() >= 0);
            best_right[read_id] = (active_right_reads.front() - right_reads.begin());
            pthread_mutex_unlock(&mutex_best_right);
        }
        // did not find valid pairs
        if (! found_pairs) {
            active_left_reads.clear();
            for (vector<alignment>::iterator lv_idx = left_reads.begin(); lv_idx != left_reads.end(); lv_idx++) {
                if (conf->pre_filter && (ignore_reads_left.find(lv_idx) != ignore_reads_left.end()))
                    continue;
                active_left_reads.push_back(lv_idx);
            }
            active_right_reads.clear();
            for (vector<alignment>::iterator rv_idx = right_reads.begin(); rv_idx != right_reads.end(); rv_idx++) {
                if (conf->pre_filter && (ignore_reads_right.find(rv_idx) != ignore_reads_right.end()))
                    continue;
                active_right_reads.push_back(rv_idx);
            }
        }
    } else {
        for (vector<alignment>::iterator lv_idx = left_reads.begin(); lv_idx != left_reads.end(); lv_idx++) {
            if (conf->pre_filter && (ignore_reads_left.find(lv_idx) != ignore_reads_left.end()))
                continue;
            active_left_reads.push_back(lv_idx);
        }
        for (vector<alignment>::iterator rv_idx = right_reads.begin(); rv_idx != right_reads.end(); rv_idx++) {
            if (conf->pre_filter && (ignore_reads_right.find(rv_idx) != ignore_reads_right.end()))
                continue;
            active_right_reads.push_back(rv_idx);
        }
    }
}

void process_data_online(ThreadData_Online* data) {

    vector<vector<alignment>::iterator> active_left_reads;
    vector<vector<alignment>::iterator> active_right_reads;

    vector<vector<alignment>::iterator>::iterator lv_idx;
    vector<vector<alignment>::iterator>::iterator rv_idx;
    vector<vector<alignment>::iterator>::iterator curr_best;

    vector<alignment>::iterator best_left_idx;
    vector<alignment>::iterator best_right_idx;

    set<vector<alignment>::iterator> ignore_idx_left;
    set<vector<alignment>::iterator> ignore_idx_right;

    unsigned int num_changed = 0;
    bool found_pairs = false;

    if (conf->pre_filter) {
        ignore_idx_left = filter_alignments(data->left_reads);
        pthread_mutex_lock(&mutex_best_left);
        if (ignore_idx_left.find(data->left_reads.begin() + genData->best_left[data->last_id]) != ignore_idx_left.end()) {
            (data->left_reads.begin() + genData->best_left[data->last_id])->is_best = false;
            for (vector<alignment>::iterator v_idx = data->left_reads.begin(); v_idx != data->left_reads.end(); v_idx++) {
                if (ignore_idx_left.find(v_idx) == ignore_idx_left.end()) {
                    genData->best_left[data->last_id] = v_idx - data->left_reads.begin();
                    v_idx->is_best = true;
                    break;
                }
            }
        }
        pthread_mutex_unlock(&mutex_best_left);

        ignore_idx_right = filter_alignments(data->right_reads);
        pthread_mutex_lock(&mutex_best_right);
        if (ignore_idx_right.find(data->right_reads.begin() + genData->best_right[data->last_id]) != ignore_idx_right.end()) {
            (data->right_reads.begin() + genData->best_right[data->last_id])->is_best = false;
            for (vector<alignment>::iterator v_idx = data->right_reads.begin(); v_idx != data->right_reads.end(); v_idx++) {
                if (ignore_idx_right.find(v_idx) == ignore_idx_right.end()) {
                    genData->best_right[data->last_id] = v_idx - data->right_reads.begin();
                    v_idx->is_best = true;
                    break;
                }
            }
        }
        pthread_mutex_unlock(&mutex_best_right);
    }

    get_active_reads(data->last_id, ignore_idx_left, ignore_idx_right, data->left_reads, data->right_reads, active_left_reads, active_right_reads, genData->best_left, genData->best_right, found_pairs);

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
        if (! best_found) 
            assert(best_found);
        bool changed = false;
        for(lv_idx = active_left_reads.begin(), rv_idx = active_right_reads.begin(); lv_idx < active_left_reads.end() && rv_idx  < active_right_reads.end(); lv_idx++, rv_idx++) {
            if ((*lv_idx)  == best_left_idx && (*rv_idx) == best_right_idx)
                continue;
            if (compare_pair(*lv_idx, *rv_idx, best_left_idx, best_right_idx, genData->coverage_map, conf)) {
                changed = true;
                best_left_idx->is_best = false;
                best_right_idx->is_best = false;
                update_coverage_map(best_left_idx, genData->coverage_map[best_left_idx->chr], 0);
                update_coverage_map(best_right_idx, genData->coverage_map[best_left_idx->chr], 0);
                best_left_idx = *lv_idx;
                best_right_idx = *rv_idx;
                best_left_idx->is_best = true;
                best_right_idx->is_best = true;
                update_coverage_map(best_left_idx, genData->coverage_map[best_left_idx->chr], 1);
                update_coverage_map(best_right_idx, genData->coverage_map[best_left_idx->chr], 1);
                
                pthread_mutex_lock(&mutex_best_left);
                genData->best_left[data->last_id] = (best_left_idx - data->left_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_left);
                pthread_mutex_lock(&mutex_best_right);
                genData->best_right[data->last_id] = (best_right_idx - data->right_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_right);
            }
        }
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
        if (active_left_reads.size() > 0)
            assert(best_found);
        bool changed = false;
        for (lv_idx = active_left_reads.begin(); lv_idx != active_left_reads.end(); lv_idx++) {
            if (lv_idx == curr_best)
                continue;

            // check if lv_idx < curr_best
            if (compare_single(*lv_idx, *curr_best, genData->coverage_map, conf)) {
                changed = true;
                (*curr_best)->is_best = false;
                update_coverage_map(*curr_best, genData->coverage_map[(*curr_best)->chr], 0);
                curr_best = lv_idx;          
                (*curr_best)->is_best = true;
                update_coverage_map(*curr_best, genData->coverage_map[(*curr_best)->chr], 1);
                pthread_mutex_lock(&mutex_best_left);
                genData->best_left[data->last_id] = (*curr_best - data->left_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_left);
            }
        }
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
        if (active_right_reads.size() > 0)
            assert(best_found);
        changed = false;
        for (rv_idx = active_right_reads.begin(); rv_idx != active_right_reads.end(); rv_idx++) {
            if (rv_idx == curr_best)
                continue;

            // check if rv_idx < curr_best
            if (compare_single(*rv_idx, *curr_best, genData->coverage_map, conf)) {
                changed = true;
                (*curr_best)->is_best = false;
                update_coverage_map(*curr_best, genData->coverage_map[(*curr_best)->chr], 0);
                curr_best = rv_idx;          
                (*curr_best)->is_best = true;
                update_coverage_map(*curr_best, genData->coverage_map[(*curr_best)->chr], 1);
                pthread_mutex_lock(&mutex_best_right);
                genData->best_right[data->last_id] = (*curr_best - data->right_reads.begin()); 
                pthread_mutex_unlock(&mutex_best_right);
            }
        }
        if (changed) 
            num_changed++;
    }
    
    delete data;

}

void *process_data_online_wrapper(void *arg) {

    while (true) {
        pthread_mutex_lock(&mutex_done);
        pthread_mutex_lock (&mutex_fifo);
        if (done && fifo.size() == 0) {
            pthread_mutex_unlock(&mutex_done);
            pthread_mutex_unlock (&mutex_fifo);
            break;
        }

        if (fifo.size() > 0) {
            ThreadData_Online* data;
            data = fifo.front();
            fifo.pop();
            pthread_mutex_unlock (&mutex_fifo);
            pthread_mutex_unlock (&mutex_done);
            pthread_cond_signal (&fifo_is_full);
            process_data_online(data);
        } else {
            pthread_mutex_unlock (&mutex_fifo);
            pthread_mutex_unlock (&mutex_done);
        }
    }
    pthread_exit((void*) 0);

}

char* parse_file_by_read(FILE* infile, ThreadData_Online* data, char* last_line, unsigned int &counter) {

    char line[1000] ;
    char cp_line[1000];
    char* ret = last_line;

    unsigned char pair_info = 0;

    alignment curr_alignment;
    string id;
    data->last_id.clear();

    unordered_map <string, size_t, hash<string> >::iterator b_idx;

    data->left_reads.clear();
    data->right_reads.clear();

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

        id = get_alignment(sl, curr_alignment, genData->chr_num, pair_info);

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
            if ((! id.compare(data->last_id)) || data->last_id.size() == 0) {
                b_idx = genData->best_left.find(id);
                if (b_idx == genData->best_left.end()) {
                    genData->best_left.insert(pair<string, size_t>(id, data->left_reads.size()));
                    curr_alignment.is_best = true;
                }
                else if (b_idx->second == data->left_reads.size())
                    curr_alignment.is_best = true;
                data->left_reads.push_back(curr_alignment);
                data->last_id = id;
            } else {
                break ;
            }
        } else {
            // id == last_id or last_id is empty
            if ((! id.compare(data->last_id)) || data->last_id.size() == 0) {
                b_idx = genData->best_right.find(id);
                if (b_idx == genData->best_right.end()) {
                    genData->best_right.insert(pair<string, size_t>(id, data->right_reads.size()));
                    curr_alignment.is_best = true;
                }
                else if (b_idx->second == data->right_reads.size())
                    curr_alignment.is_best = true;
                data->right_reads.push_back(curr_alignment);
                data->last_id = id;
            } else {
                break ;
            }
        }
    }
    strcpy(last_line, cp_line);
    return ret;
}

void parse_complete_file(ThreadData_Batch* data) {

    char line[1000] ;
    char cp_line[1000];
    char* ret;

    unsigned int counter = 0;
    unsigned char pair_info = 0;
    
    alignment curr_alignment;
    string id;
    string last_id = string("");

    FILE* infile = open_bam_pipe_in(conf->infile);
    bool read_header = false;

    if (conf->verbose) {
        fprintf(stdout, "\nReading input file from: %s\n", conf->infile.c_str());
    }

    while (true) {

        ret = fgets(line, sizeof(line), infile);
        strcpy(cp_line, line);
        counter++;

        if (!ret) {
            if (counter == 1) {
                fprintf(stderr, "Could not read SAM file %s\n", conf->infile.c_str());
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

        id = get_alignment(sl, curr_alignment, genData->chr_num, pair_info);

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
            unordered_map<string, vector<alignment> >::iterator id_idx = data->read_map_left.find(id);
            if (id_idx != data->read_map_left.end()) {
                id_idx->second.push_back(curr_alignment);
            } else {
                // create new element in read map
                vector<alignment> tmp_align;
                tmp_align.push_back(curr_alignment);
                data->read_map_left.insert( pair<string, vector<alignment> >(id, tmp_align) );
            }
        } else {
            unordered_map<string, vector<alignment> >::iterator id_idx = data->read_map_right.find(id);
            if (id_idx != data->read_map_right.end()) {
                id_idx->second.push_back(curr_alignment);
            } else {
                // create new element in read map
                vector<alignment> tmp_align;
                tmp_align.push_back(curr_alignment);
                data->read_map_right.insert( pair<string, vector<alignment> >(id, tmp_align) );
            }
        }
    }

    if (conf->verbose) 
        fprintf(stdout, "\nsuccessfully parsed %i lines\n", counter);

    fclose(infile);
}

int main(int argc, char *argv[]) {

    conf = new Config::Config(argc, argv);

    if (conf->parse_complete) {
        ThreadData_Batch* data = new ThreadData_Batch::ThreadData_Batch();
        parse_complete_file(data);

        // pre filter the alignments
        if (conf->pre_filter) {
            pre_filter_alignment_map(data->read_map_left);
            pre_filter_alignment_map(data->read_map_right);
        }

        if (conf->verbose) 
            fprintf(stdout, "\nBuilding up the coverage map ...\n");

        // build coverage map
        unordered_map<string, vector<alignment> >::iterator r_idx;
        unordered_map<string, vector<alignment> >::iterator l_idx;
        vector<alignment>::iterator curr_best;
        unsigned char max_qual;
        
        // handle left reads
        for (l_idx = data->read_map_left.begin(); l_idx != data->read_map_left.end(); l_idx++) {
            max_qual = 0;
            for (vector<alignment>::iterator v_idx = l_idx->second.begin(); v_idx != l_idx->second.end(); v_idx++) {
                if (v_idx->quality > max_qual) {
                    curr_best = v_idx;
                    max_qual = curr_best->quality;
                }
            }
            if (max_qual == 0) 
                curr_best = l_idx->second.begin();
            curr_best->is_best = true;
            update_coverage_map(curr_best, genData->coverage_map[curr_best->chr], 1);
        }
        // handle right reads
        for (r_idx = data->read_map_right.begin(); r_idx != data->read_map_right.end(); r_idx++) {
            max_qual = 0;
            for (vector<alignment>::iterator v_idx = r_idx->second.begin(); v_idx != r_idx->second.end(); v_idx++) {
                if (v_idx->quality > max_qual) {
                    curr_best = v_idx;
                    max_qual = curr_best->quality;
                }
            }
            if (max_qual == 0) 
                curr_best = r_idx->second.begin();
            curr_best->is_best = true;
            update_coverage_map(curr_best, genData->coverage_map[curr_best->chr], 1);
        }
        if (conf->verbose) 
            fprintf(stdout, "... done.\n");

        if (conf->verbose) 
            fprintf(stdout, "\nSmoothing the coverage map ...\n");

        double sum_min_loss = 0.0;
        unsigned int num_ambiguous_single = 0;
        unsigned int num_ambiguous_paired = 0;
        int num_changed_single = 0;
        int num_changed_paired = 0;

        // computing active set of reads
        get_active_read_set(data);

        if (conf->verbose) { 
            sum_min_loss = get_total_min_loss(data->read_map_left, genData->coverage_map);
            sum_min_loss += get_total_min_loss(data->read_map_right, genData->coverage_map);
            fprintf(stdout, "\nObjective before smoothing:\n\tchanged single (paired):  %i/%i (%i/%i)\n\tobjective: %lf\n", num_changed_single, num_ambiguous_single, num_changed_paired, num_ambiguous_paired, sum_min_loss); 
        }

        // iterate to smooth coverage map
        for (unsigned int iteration = 0; iteration < conf->iterations; iteration++) {

            if (conf->verbose) 
                fprintf(stdout, "\n\t... processing iteration %i of %i ...\n", iteration + 1, conf->iterations);

            num_ambiguous_single = 0;
            num_ambiguous_paired = 0;

            if (conf->use_pair_info) {
                num_changed_paired = smooth_coverage_map_paired(data->active_left_pair, data->active_right_pair, genData->coverage_map, num_ambiguous_paired);
            } 
            num_changed_single = smooth_coverage_map_single(data->active_left_single, genData->coverage_map, num_ambiguous_single);
            num_changed_single += smooth_coverage_map_single(data->active_right_single, genData->coverage_map, num_ambiguous_single);

            if (conf->verbose) { 
                sum_min_loss = get_total_min_loss(data->read_map_left, genData->coverage_map);
                sum_min_loss += get_total_min_loss(data->read_map_right, genData->coverage_map);
                fprintf(stdout, "\n\tchanged single (paired):  %i/%i (%i/%i)\n\tobjective: %lf\n", num_changed_single, num_ambiguous_single, num_changed_paired, num_ambiguous_paired, sum_min_loss); 
            }

            if ((num_changed_single + num_changed_paired) == 0) {
                if (conf->verbose)
                    fprintf(stdout, "\n\tNo further improvement expected - leave iterations now.\n"); 
                break;
            }
        }

        if (conf->verbose)  
            fprintf(stdout, "... done.\n\n");

        write_output(data->read_map_left, data->read_map_right, genData->chr_num);
    } else {
        
        pthread_t* threads = NULL;
        unsigned int tid;

        for (unsigned int iteration = 0; iteration < conf->iterations; iteration++) {
            
            ThreadData_Online* data = new ThreadData_Online::ThreadData_Online();
            pthread_mutex_lock(&mutex_done);
            done = false;
            pthread_mutex_unlock(&mutex_done);

            // check for multi-threading
            if (conf->num_threads > 1) {
                // create threads
                threads = new pthread_t[conf->num_threads - 1];
                
                pthread_attr_init(&attr);
                pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                
                for (tid = 0; tid < conf->num_threads - 1; tid++) {
                   pthread_create(&threads[tid], &attr, process_data_online_wrapper, (void *) tid); 
                }
            }

            char line[1000];
            char last_line[1000];

            FILE* infile = open_bam_pipe_in(conf->infile);

            if (! infile) {
                fprintf(stderr, "Could not open %s for reading!\n", conf->infile.c_str());
                exit(1);
            }
            char* ret = fgets(line, 1000, infile);
            unsigned int counter = 0;

            if (!ret) {
                fprintf(stderr, "Could not read SAM file %s\n", conf->infile.c_str());
                exit(1);
            }

            if (conf->verbose)
                fprintf(stdout, "\nReading input file from: %s (iteration %i)\n", conf->infile.c_str(), iteration + 1);

            if (iteration == 0) {
                
                // search for header
                char* sl = strtok(line, "\t");
                bool header_parsed = false;
                while (line[0] == '@') {
                    parse_header(sl);
                    ret = fgets(line, sizeof(line), infile);
                    strcpy(last_line, line);
                    sl = strtok(line, "\t");
                    header_parsed = true;
                }

                if (! header_parsed) {
                    fprintf(stderr, "No header found in %s.\n", conf->infile.c_str());
                    exit(1);
                }
            
            } else {
                while (line[0] == '@') {
                    ret = fgets(line, sizeof(line), infile);
                }
                strcpy(last_line, line);
            }

            string last_id = string("");
            // fill FIFO
            while ((ret = parse_file_by_read(infile, data, last_line, counter)) || data->left_reads.size() > 0 || data->right_reads.size() > 0) {

                pthread_mutex_lock (&mutex_fifo);
                if (fifo.size() >= conf->max_fifo_size) pthread_cond_wait(&fifo_is_full, &mutex_fifo);
                pthread_mutex_unlock (&mutex_fifo);
                        
                if (conf->num_threads < 2) {
                    process_data_online(data);
                } else {
                    pthread_mutex_lock (&mutex_fifo);
                    fifo.push(data);
                    pthread_mutex_unlock (&mutex_fifo);
                }
                
                if (!ret)
                    break;
                else
                    data = new ThreadData_Online::ThreadData_Online();
            }

            fclose(infile);
            pthread_mutex_lock(&mutex_done);
            done = true;
            pthread_mutex_unlock(&mutex_done);
            // join threads, if neccessary
            if (conf->num_threads > 1) {
                fprintf(stdout, "waiting for threads to join\n");
                for (tid = 0; tid < conf->num_threads - 1; tid++) {
                   pthread_join(threads[tid], NULL); 
                }
            }

            if (conf->verbose) 
                fprintf(stdout, "\nsuccessfully parsed %i lines\n", counter - 1);

        }

        write_output_direct();

        if (conf->num_threads > 1) delete threads;
    }

    delete conf;
    delete genData;
   

    return 0;
}

