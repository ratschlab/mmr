#include <string>
#include <cstdio>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>

#include "Alignment.h"
#include "GeneralData.h"
#include "config.h"

using namespace std;

extern GeneralData* genData;
extern Config* conf;

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

double get_variance(vector<unsigned short> &exon_coverage, vector<unsigned short> &intron_coverage) {

    if (exon_coverage.size() < 2) {
        return 0.0;
    } else {
        double sum = 0.0;
        double int_pnlty = 0.0;
        for ( vector<unsigned short>::iterator it = exon_coverage.begin(); it != exon_coverage.end(); it++) {
            sum += (double) *it;
        }
        double mean = sum / (double) exon_coverage.size();
        sum = 0.0;
        for ( vector<unsigned short>::iterator it = exon_coverage.begin(); it != exon_coverage.end(); it++) {
            sum += pow((double) *it - mean, 2.0);
        }

        if (intron_coverage.size() > 0) {
            int_pnlty = intron_penalty(intron_coverage);
        }

        return (sum / ((double) exon_coverage.size() - 1.0)) + int_pnlty;  
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

bool compare_pair(vector<Alignment>::iterator candidate_left, vector<Alignment>::iterator candidate_right, vector<Alignment>::iterator best_left, vector<Alignment>::iterator best_right) {

    bool used_mip = false;
    double candidate_loss = 0.0;
    double best_loss = 0.0;

    // check how objective is to be computed
    // we need segments overlapping at least one of the two candidates
    if (conf->use_mip_objective) {
       
        // compute overlap between candidate alignments
        set<unsigned long> overlap;
        set<unsigned long> overlap_cl_bl = candidate_left->get_overlap(*best_left);
        set<unsigned long> overlap_cl_br = candidate_left->get_overlap(*best_right); 
        set<unsigned long> overlap_cr_bl = candidate_right->get_overlap(*best_left);
        set<unsigned long> overlap_cr_br = candidate_right->get_overlap(*best_right); 

        // loss candidate_left
        set_union(overlap_cl_bl.begin(), overlap_cl_bl.end(), overlap_cl_br.begin(), overlap_cl_br.end(), inserter(overlap, overlap.begin()));
        pair<double, double> loss_cl = genData->segments.get_exon_segment_loss(candidate_left, overlap, false);
        overlap.clear();

        // loss candidate_rigth
        set_union(overlap_cr_bl.begin(), overlap_cr_bl.end(), overlap_cr_br.begin(), overlap_cr_br.end(), inserter(overlap, overlap.begin()));
        pair<double, double> loss_cr = genData->segments.get_exon_segment_loss(candidate_right, overlap, false);
        overlap.clear();

        // loss best_left
        set_union(overlap_cr_bl.begin(), overlap_cr_bl.end(), overlap_cl_bl.begin(), overlap_cl_bl.end(), inserter(overlap, overlap.begin()));
        pair<double, double> loss_bl = genData->segments.get_exon_segment_loss(best_left, overlap, true);
        overlap.clear();

        // loss best_right
        set_union(overlap_cr_br.begin(), overlap_cr_br.end(), overlap_cl_br.begin(), overlap_cl_br.end(), inserter(overlap, overlap.begin()));
        pair<double, double> loss_br = genData->segments.get_exon_segment_loss(best_right, overlap, true);
        overlap.clear();

        // check if any loss is valid
        // first is loss_with and second is loss_without
        if (loss_cl.first >= 0.0 || loss_cr.first >= 0.0 || loss_br.first >= 0.0 || loss_bl.first >= 0.0) {
            candidate_loss +=  (loss_cl.first >= 0.0) ? loss_cl.first : 0;
            candidate_loss +=  (loss_cr.first >= 0.0) ? loss_cl.first : 0;
            candidate_loss +=  (loss_bl.second >= 0.0) ? loss_bl.second : 0;
            candidate_loss +=  (loss_br.second >= 0.0) ? loss_br.second : 0;
            best_loss +=  (loss_cl.second >= 0.0) ? loss_cl.second : 0;
            best_loss +=  (loss_cr.second >= 0.0) ? loss_cl.second : 0;
            best_loss +=  (loss_bl.first >= 0.0) ? loss_bl.first : 0;
            best_loss +=  (loss_br.first >= 0.0) ? loss_br.first : 0;
            used_mip = true;
        }
    }

    if (! conf->use_mip_objective || ! used_mip) {
        vector<unsigned short> intron_cov_candidate_left;
        vector<unsigned short> intron_cov_candidate_right;
        vector<unsigned short> intron_cov_best_left;
        vector<unsigned short> intron_cov_best_right;

        vector<unsigned short> exon_cov_candidate_left_without;
        candidate_left->get_coverage(conf->window_size, exon_cov_candidate_left_without, intron_cov_candidate_left);
        vector<unsigned short> exon_cov_candidate_left_with = alter_coverage(exon_cov_candidate_left_without, min((unsigned int) candidate_left->start, conf->window_size), conf->window_size, true);

        vector<unsigned short> exon_cov_candidate_right_without;
        candidate_right->get_coverage(conf->window_size, exon_cov_candidate_right_without, intron_cov_candidate_right);
        vector<unsigned short> exon_cov_candidate_right_with =alter_coverage(exon_cov_candidate_right_without, min((unsigned int) candidate_right->start, conf->window_size), conf->window_size, true);

        vector<unsigned short> exon_cov_best_left_with;
        best_left->get_coverage(conf->window_size, exon_cov_best_left_with, intron_cov_best_left);
        vector<unsigned short> exon_cov_best_left_without = alter_coverage(exon_cov_best_left_with, min((unsigned int) best_left->start, conf->window_size), conf->window_size, false);
        vector<unsigned short> exon_cov_best_right_with;
        best_right->get_coverage(conf->window_size, exon_cov_best_right_with, intron_cov_best_right);
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

        best_loss = var_cand_left_without + var_cand_right_without + var_best_left_with + var_best_right_with;
        candidate_loss = var_cand_left_with + var_cand_right_with + var_best_left_without + var_best_right_without;
    }
    return candidate_loss < best_loss;
}

bool compare_pair_mip(vector<Alignment>::iterator candidate_left, vector<Alignment>::iterator candidate_right, vector<Alignment>::iterator best_left, vector<Alignment>::iterator best_right) {

    vector<unsigned short> intron_cov_candidate_left;
    vector<unsigned short> intron_cov_candidate_right;
    vector<unsigned short> intron_cov_best_left;
    vector<unsigned short> intron_cov_best_right;

    vector<unsigned short> exon_cov_candidate_left_without;
    candidate_left->get_coverage(conf->window_size, exon_cov_candidate_left_without, intron_cov_candidate_left);
    vector<unsigned short> exon_cov_candidate_left_with = alter_coverage(exon_cov_candidate_left_without, min((unsigned int) candidate_left->start, conf->window_size), conf->window_size, true);

    vector<unsigned short> exon_cov_candidate_right_without;
    candidate_right->get_coverage(conf->window_size, exon_cov_candidate_right_without, intron_cov_candidate_right);
    vector<unsigned short> exon_cov_candidate_right_with =alter_coverage(exon_cov_candidate_right_without, min((unsigned int) candidate_right->start, conf->window_size), conf->window_size, true);

    vector<unsigned short> exon_cov_best_left_with;
    best_left->get_coverage(conf->window_size, exon_cov_best_left_with, intron_cov_best_left);
    vector<unsigned short> exon_cov_best_left_without = alter_coverage(exon_cov_best_left_with, min((unsigned int) best_left->start, conf->window_size), conf->window_size, false);
    vector<unsigned short> exon_cov_best_right_with;
    best_right->get_coverage(conf->window_size, exon_cov_best_right_with, intron_cov_best_right);
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


bool compare_single(vector<Alignment>::iterator candidate, vector<Alignment>::iterator best) {

    vector<unsigned short> intron_cov_candidate;
    vector<unsigned short> intron_cov_best;
    vector<unsigned short> exon_cov_candidate_without; 
    candidate->get_coverage(conf->window_size, exon_cov_candidate_without, intron_cov_candidate);
    vector<unsigned short> exon_cov_best_with;
    best->get_coverage(conf->window_size, exon_cov_best_with, intron_cov_best);
    vector<unsigned short> exon_cov_candidate_with = alter_coverage(exon_cov_candidate_without, min((unsigned int) candidate->start, conf->window_size), conf->window_size, true);
    vector<unsigned short> exon_cov_best_without = alter_coverage(exon_cov_best_with, min((unsigned int) best->start, conf->window_size), conf->window_size, false);

    vector<unsigned short> empty_cov;
    double var_cand_without = get_variance(exon_cov_candidate_without, empty_cov);
    double var_best_with = get_variance(exon_cov_best_with, intron_cov_best);
    double var_cand_with = get_variance(exon_cov_candidate_with, intron_cov_candidate);
    double var_best_without = get_variance(exon_cov_best_without, empty_cov);

    return (var_cand_with + var_best_without) < (var_cand_without + var_best_with);
}

set<vector<Alignment>::iterator> filter_alignments(vector<Alignment> &aligns) {

    set<vector<Alignment>::iterator> to_erase;
    vector<Alignment>::iterator v_idx;
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

void get_plifs_from_file(string &filename) {
    
    FILE* infile = fopen(filename.c_str(), "r");    
    char* ret;
    char line[1000];
    int i_count;
    double c, sl1, sl2, sr1, sr2;

    map< double, vector<double> > plifs;

    while (true) {
        ret = fgets(line, sizeof(line), infile);

        if (!ret)
            break;
        
        if (line[0] == '#')
            continue;

        i_count = sscanf(line, "%f %f %f %f %f", &c, &sl1, &sl2, &sr1, &sr2);

        if (i_count != 5) {
            fprintf(stderr, "Number of items per line in PLIF-Input file is different from 5! - Bail out!\n");
            exit(1);
        }

        vector<double> tmp_plif;
        tmp_plif.push_back(sl1);
        tmp_plif.push_back(sl2);
        tmp_plif.push_back(sr1);
        tmp_plif.push_back(sr2);
        plifs.insert(plifs.begin(), pair<double, vector<double> >(c, tmp_plif));
    }

    genData->plifs = plifs;
    
}

double compute_loss(double observed_cov, double predicted_cov) {
    
    // get plif iterators to interpolate entries
    map< double, vector<double> >::iterator lower = genData->plifs.upper_bound(observed_cov);
    lower--;
    map< double, vector<double> >::iterator upper = genData->plifs.upper_bound(observed_cov);

    double loss = 0.0;

    if (observed_cov < 0)
        observed_cov = 0;

    // determine if left or right part of function needs to be queried
    if (predicted_cov > observed_cov) {
        if (observed_cov > 30000)
            observed_cov = 30000;
        loss = ((lower->second.at(1) + upper->second.at(1)) / 2)*predicted_cov*predicted_cov + ((lower->second.at(2) + upper->second.at(2)) / 2)*predicted_cov; 
    } else {
        if (observed_cov > 30000)
            observed_cov = 30000;
        loss = ((lower->second.at(3) + upper->second.at(3)) / 2)*predicted_cov*predicted_cov + ((lower->second.at(4) + upper->second.at(4)) / 2)*predicted_cov; 
    }

    return loss;
}

