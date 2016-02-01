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

#include <string>
#include <cstdio>
#include <vector>
#include <set>
#include <deque>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <unistd.h>

#include "Alignment.h"
#include "GeneralData.h"
#include "config.h"
#include "Utils.h" 
#include "bam_sort.h"

using namespace std;

extern GeneralData* genData;
extern Config* conf;

FILE* open_bam_pipe_in(std::string & in_fname)
{
	string command = conf->samtools + string(" view -h 2> /dev/null ") + in_fname  + " && echo samtools subprocess for reading terminated successfully";
	FILE* IN_FP=NULL ;
	fflush(stdout) ;
    for (int i = 0; i < 10; i++) {
        IN_FP = popen(command.c_str(), "r") ;
        if (IN_FP)
            break;
        sleep(1);
    }
	return IN_FP ;
}

FILE* open_bam_pipe_out(std::string & out_fname)
{
	string command = conf->samtools + string(" view -bS -o ") + out_fname +  " - " + " && echo samtools for writing subprocess terminated successfully";
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

double intron_penalty(vector<unsigned int> &intron_coverage) {
    double sum = 0.0;

    if (intron_coverage.size() > 0) {
        for (vector<unsigned int>::iterator it = intron_coverage.begin(); it != intron_coverage.end(); it++) {
            sum += (double) *it;
        }
        sum /= (double) intron_coverage.size();
    }

    return sum;
}


double get_variance_global() {

    double total_var = 0.0;
    double total_cov = 0.0;

    // iterate over chromosomes
    for (map <pair<unsigned int, unsigned char>, vector<unsigned int> >::iterator it = genData->coverage_map.begin(); it != genData->coverage_map.end(); it++) {

        double sum = 0.0;
        double sum_sq = 0.0;
        double mean = 0.0;
        double var = 0.0;
        deque<double> cov;
        
        // collect segment sizes in case we have annotated segment information
        vector<unsigned long> seg_size;
        if (conf->use_brkpts) {
            unsigned long s = 1;
            for (size_t i = 1; i < genData->breakpoint_map[it->first].size(); i++) {
                if (genData->breakpoint_map[it->first].at(i)) {
                    seg_size.push_back(s);
                    s = 0;
                }
                s++;
            }
            if (s > 0)
                seg_size.push_back(s);
        } else {
            seg_size.push_back(it->second.size());
        }   

        size_t cov_idx = 0;
        for (size_t j = 0; j < seg_size.size(); j++) {
            if (seg_size.at(j) < 2) {
                cov_idx++;
                continue;
            }
            sum = 0.0;
            sum_sq = 0.0;
            mean = 0.0;
            var = 0.0;
            cov.clear();
            unsigned long max_size = min((unsigned long) conf->window_size, seg_size.at(j));
            for (size_t i = 0; i < seg_size.at(j); i++) {
                total_cov += it->second.at(cov_idx);
                if (i < max_size) {
                    cov.push_back((double) it->second.at(cov_idx));
                    sum += (double) it->second.at(cov_idx);
                    sum_sq += pow((double) it->second.at(cov_idx), 2.0);
                    cov_idx++;
                    if (i == (max_size - 1)) {
                        mean = sum / (double) cov.size();
                        var = (sum_sq + (mean * (double) cov.size() * mean) - (mean * 2 * sum));
                    } else {
                        continue;
                    }
                } else {
                    cov.push_back((double) it->second.at(cov_idx));
                    cov_idx++;
                    sum += (cov.back() - cov.front());
                    mean = sum / ((double) cov.size() - 1.0);
                    sum_sq += (pow(cov.back(), 2.0) - pow(cov.front(), 2.0));
                    cov.pop_front();
                    var = (sum_sq + (mean * (double) cov.size() * mean) - (mean * 2 * sum));
                }
                total_var += (var / ((double) cov.size() - 1.0));
                //total_var += (sqrt(var / ((double) cov.size() - 1.0)) / mean);
                //if (it->first.first == 5 && cov_idx >= 107552486 && cov_idx <= 107556023) { 
                //    fprintf(stdout, "glob (%f): ", var / ((double) cov.size() - 1.0));
                //    for (deque<double>::iterator c = cov.begin(); c != cov.end(); c++)
                //        fprintf(stdout, "%.0f ", *c);
                //    fprintf(stdout, "\n");
                //}
            }
        }
    }
    //fprintf(stdout, "total cov %f\n", total_cov);
    return total_var;
}

double get_variance(vector<vector<vector<unsigned long> > > &exon_coverage, vector<vector<set<unsigned long> > > &genome_pos) {

    double total_var = 0.0;
    set<unsigned long> already_seen;
    set<unsigned long>::iterator gp;
    for (size_t k = 0; k < exon_coverage.size(); k++) {
        for (size_t j = 0; j < exon_coverage.at(k).size(); j++) {
            if (exon_coverage.at(k).at(j).size() < 2)
                continue;
            double sum = 0.0;
            double sum_sq = 0.0;
            double mean = 0.0;
            double var = 0.0;
            deque<double> cov;
            unsigned long max_size = min((unsigned long) conf->window_size, exon_coverage.at(k).at(j).size());
            for (size_t i = 0; i < exon_coverage.at(k).at(j).size(); i++) {
                if (i == 0)
                    gp = genome_pos.at(k).at(j).begin();
                else
                    gp++;
                if (i < max_size) {
                    cov.push_back((double) exon_coverage.at(k).at(j).at(i));
                    sum += (double) exon_coverage.at(k).at(j).at(i);
                    sum_sq += pow((double) exon_coverage.at(k).at(j).at(i), 2.0);
                    if (i == (max_size - 1)) {
                        mean = sum / (double) cov.size();
                        var = (sum_sq + (mean * (double) cov.size() * mean) - (mean * 2 * sum));
                    } else {
                        continue;
                    }
                } else {
                    cov.push_back((double) exon_coverage.at(k).at(j).at(i));
                    sum += (cov.back() - cov.front());
                    mean = sum / ((double) cov.size() - 1.0);
                    sum_sq += (pow(cov.back(), 2.0) - pow(cov.front(), 2.0));
                    cov.pop_front();
                    var = (sum_sq + (mean * (double) cov.size() * mean) - (mean * 2 * sum));
                }
                if (!(already_seen.find(*gp) != already_seen.end())) {
                    already_seen.insert(*gp);
                    total_var += (var / ((double) cov.size() - 1.0));
                }
                //fprintf(stdout, "loc (%f): ", var / ((double) cov.size() - 1.0));
                //for (deque<double>::iterator c = cov.begin(); c != cov.end(); c++)
                //    fprintf(stdout, "%.0f ", *c);
                //fprintf(stdout, "\n");
            }
        }
    }
    return total_var;
}


vector<unsigned int> alter_coverage(vector<unsigned int> &source, unsigned int window_left, unsigned int window_right, bool is_positive) {

    vector<unsigned int> target;

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

    char cp_line[10000];
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
    delete sl;
    return return_line.substr(1, return_line.size() - 2);
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
                genData->chr_num.insert( pair<string, unsigned int>(tmp_sl,  (unsigned int) genData->chr_num.size() + 1) );
            }
            else {
                fprintf(stderr, "WARNING: Doubled contig names in header!\n Ignoring %s\n\n", tmp_sl.c_str());
            }
        }
        else if (idx == 2) {
            if (genData->chr_size_cum.size() == 0)
                genData->chr_size_cum.push_back(0);
            else
                genData->chr_size_cum.push_back(genData->chr_size_cum.back() + genData->chr_size.back());
            tmp_sl = tmp_sl.substr(3, tmp_sl.size());
            genData->chr_size.push_back(atoi(tmp_sl.c_str()));    
            vector<unsigned int> tmp_cov(genData->chr_size.back(), 0);
            if (conf->verbose) 
                fprintf(stdout, "\t...reserving memory for contig %s (+) of size %i\n", chr_name.c_str(), genData->chr_size.back());
            genData->coverage_map.insert( pair<pair<unsigned int, unsigned char>, vector<unsigned int> >(pair<unsigned int, unsigned char>(genData->chr_num[chr_name], '+'), tmp_cov) ); 
            if (conf->use_brkpts) {
                vector<bool> tmp_brkp(genData->chr_size.back(), 0);
                genData->breakpoint_map.insert( pair<pair<unsigned int, unsigned char>, vector<bool> >(pair<unsigned int, unsigned char>(genData->chr_num[chr_name], '+'), tmp_brkp) ); 
            } else {
                vector<bool> tmp_brkp;
                genData->breakpoint_map.insert( pair<pair<unsigned int, unsigned char>, vector<bool> >(pair<unsigned int, unsigned char>(genData->chr_num[chr_name], '+'), tmp_brkp) ); 
            }
            if (conf->strand_specific) {
                if (conf->verbose) 
                    fprintf(stdout, "\t...reserving memory for contig %s (-) of size %i\n", chr_name.c_str(), genData->chr_size.back());
                genData->coverage_map.insert( pair<pair<unsigned int, unsigned char>, vector<unsigned int> >(pair<unsigned int, unsigned char>(genData->chr_num[chr_name], '-'), tmp_cov) ); 
            }
        }   
        idx ++;
        sl = strtok(NULL, "\t");    
    }

    if (genData->chr_size.size() != genData->chr_num.size() || (!conf->strand_specific && genData->chr_num.size() != genData->coverage_map.size()) || (conf->strand_specific && 2*genData->chr_num.size() != genData->coverage_map.size())) {
        fprintf(stderr, "\nERROR: Header information incomplete!");
        exit(-1);
    }
}

bool compare_pair(vector<Alignment>::iterator candidate_left, vector<Alignment>::iterator candidate_right, vector<Alignment>::iterator best_left, vector<Alignment>::iterator best_right, double &loss, double &gain, bool debug) {

    bool used_mip = false;
    double candidate_loss = 0.0;
    double best_loss = 0.0;

    // check how objective is to be computed
    // we need segments overlapping at least one of the two candidates
    if (conf->use_mip_objective) {
       
        // compute overlap between candidate alignments
        set<unsigned long> overlap;
        set<unsigned long> genome_pos_cl = candidate_left->get_exon_pos();
        set<unsigned long> genome_pos_cr = candidate_right->get_exon_pos();
        set<unsigned long> genome_pos_bl = best_left->get_exon_pos();
        set<unsigned long> genome_pos_br = best_right->get_exon_pos();

        pair<double, double> loss_best, loss_cand;
        pair<double, double> loss_cl, loss_cr, loss_bl, loss_br;
        // loss candidate
        vector<vector<Alignment>::iterator> candidate;
        candidate.push_back(candidate_left);
        candidate.push_back(candidate_right);

        set_union(genome_pos_bl.begin(), genome_pos_bl.end(), genome_pos_br.begin(), genome_pos_br.end(), inserter(overlap, overlap.begin()));
        loss_cand = genData->segments.get_exon_segment_loss(candidate, overlap, debug);

        // loss best
        vector<vector<Alignment>::iterator> best;
        best.push_back(best_left);
        best.push_back(best_right);
        overlap.clear();
        set_union(genome_pos_cl.begin(), genome_pos_cl.end(), genome_pos_cr.begin(), genome_pos_cr.end(), inserter(overlap, overlap.begin()));
        loss_best = genData->segments.get_exon_segment_loss(best, overlap, debug);

        // check if any loss is valid
        // first is loss_with and second is loss_without
        if (loss_best.first >= 0.0 || loss_cand.first >= 0.0) {
           best_loss = loss_best.first + loss_cand.second;
           candidate_loss = loss_best.second + loss_cand.first;
           used_mip = true;
        }
    }

    // if we never were supposed to use the mip objective or did not use it
    // for other reasons, use the variance objective instead, but do not add it
    // to the total loss for mip-objective
    if (! conf->use_mip_objective || ! used_mip) {

        // initialize segment coverage
        vector<vector<vector<unsigned long> > > cov_keep;
        vector<vector<vector<unsigned long> > > cov_change;
        vector<vector<set<unsigned long> > > genome_pos;

        // init vectors of candidate and best alignments
        vector<pair<vector<Alignment>::iterator, bool> > aligns;
        aligns.push_back( make_pair(candidate_left, false) );
        aligns.push_back( make_pair(candidate_right, false) );
        aligns.push_back( make_pair(best_left, true) );
        aligns.push_back( make_pair(best_right, true) );

        // compute coverage
        compute_coverage_loss(aligns, cov_keep, cov_change, genome_pos);
        
        best_loss = get_variance(cov_keep, genome_pos);
        candidate_loss = get_variance(cov_change, genome_pos);

    }

    if (debug) {
        fprintf(stdout, "best loss: %f\n", best_loss);
        fprintf(stdout, "candidate loss: %f\n", candidate_loss);
        fprintf(stdout, "delta (best - cand): %f\n", best_loss - candidate_loss);
    }
 
    if (conf->use_mip_objective && ! used_mip) {
        loss = -1.0;
        gain = 0.0;
    } else {
        loss = min(candidate_loss, best_loss);
        gain = max(best_loss - candidate_loss, 0.0);
    }
    return (candidate_loss < best_loss);
}


bool compare_align_iter_start(const pair<vector<Alignment>::iterator, bool> &left, const pair<vector<Alignment>::iterator, bool> &right) {
    return left.first->start < right.first->start;
}


void compute_coverage_loss(vector<pair<vector<Alignment>::iterator,bool> > aligns, vector<vector<vector<unsigned long> > > &cov_keep, vector<vector<vector<unsigned long> > > &cov_change, vector<vector<set<unsigned long> > > &genome_pos) {

    // sort alignments by starting position
    sort(aligns.begin(), aligns.end(), compare_align_iter_start); 

    // determine genomic position set of all alignments (including windows)
    for (size_t i = 0; i < aligns.size(); i++) {
        genome_pos.push_back(aligns.at(i).first->get_genome_pos(conf->window_size));
    }

    // fill coverage maps
    for (size_t i = 0; i < aligns.size(); i++) {
        vector<vector<unsigned long> > tmp;
        cov_keep.push_back(tmp);
        aligns.at(i).first->fill_coverage_vector(cov_keep.at(i));
    }

    cov_change = cov_keep;

    // iterate over alignments and alter coverage maps
    for (size_t i = 0; i < aligns.size(); i++) {
        aligns.at(i).first->alter_coverage_vector(cov_change, genome_pos, aligns.at(i).second);
    }
}


bool compare_single(vector<Alignment>::iterator candidate, vector<Alignment>::iterator best, double &loss, double &gain, bool debug) {

    bool used_mip = false;
    double candidate_loss = 0.0;
    double best_loss = 0.0;

    pair<double, double> loss_candidate, loss_best;

    //debug=false;//true;

    // check how objective is to be computed
    // we need segments overlapping at least one of the two candidates
    if (conf->use_mip_objective) {

        // compute overlap between candidate and best alignments
        set<unsigned long> genome_pos_candidate = candidate->get_exon_pos();
        set<unsigned long> genome_pos_best = best->get_exon_pos();

        vector<vector<Alignment>::iterator> cand_tmp;
        vector<vector<Alignment>::iterator> best_tmp;
        cand_tmp.push_back(candidate);
        best_tmp.push_back(best);

        loss_candidate = genData->segments.get_exon_segment_loss(cand_tmp, genome_pos_best, debug);
        loss_best = genData->segments.get_exon_segment_loss(best_tmp, genome_pos_candidate, debug);

        // check if any loss is valid
        // first is loss_with and second is loss_without
        if (loss_candidate.first >= 0.0 || loss_best.first >= 0.0) {
            used_mip = true;
        }
        if (debug) {
            fprintf(stdout, "candidate (single):\n");
            candidate->print();
            fprintf(stdout, "best (single):\n");
            best->print();
            fprintf(stdout, "cand loss (w): %f best loss (w): %f\n", loss_candidate.first, loss_best.first);
            fprintf(stdout, "cand loss (wo): %f best loss (wo): %f\n", loss_candidate.second, loss_best.second);
        }

        candidate_loss +=  (loss_candidate.first >= 0.0) ? loss_candidate.first : 0;
        candidate_loss +=  (loss_best.second >= 0.0) ? loss_best.second : 0;
        best_loss +=  (loss_candidate.second >= 0.0) ? loss_candidate.second : 0;
        best_loss +=  (loss_best.first >= 0.0) ? loss_best.first : 0;

    }

    // if we never were supposed to use the mip objective or did not use it
    // for other reasons, use the variance objective instead, but do not add it
    // to the total loss for mip-objective
    if (! conf->use_mip_objective || ! used_mip) {

        // initialize segment coverage
        vector<vector<vector<unsigned long> > > cov_keep;
        vector<vector<vector<unsigned long> > > cov_change;
        vector<vector<set<unsigned long> > > genome_pos;

        // init vectors of candidate and best alignments
        vector<pair<vector<Alignment>::iterator, bool> > aligns;
        aligns.push_back( make_pair(candidate, false) );
        aligns.push_back( make_pair(best, true) );

        // compute coverage
        compute_coverage_loss(aligns, cov_keep, cov_change, genome_pos);
        
        best_loss = get_variance(cov_keep, genome_pos);
        candidate_loss = get_variance(cov_change, genome_pos);
        if (debug) {
            fprintf(stdout, "candidate (single):\n");
            candidate->print();
            fprintf(stdout, "best (single):\n");
            best->print();
            fprintf(stdout, "cov_keep (single):\n");
            unsigned long ss = 0;
            int pp = 0;
            for (size_t l = 0; l < cov_keep.size(); l++) {
                for (size_t k = 0; k < cov_keep.at(l).size(); k++) {
                    ss = 0;
                    pp = 0;
                    fprintf(stdout, "|");
                    for (size_t j = 0; j < cov_keep[l][k].size(); j++) {
                        ss += cov_keep[l][k][j];
                        pp++;
                        fprintf(stdout, "%lu ", cov_keep[l][k][j]);
                    }
                    //fprintf(stdout, "%i - %lu|", pp, ss);
                }
                fprintf(stdout, "\n");
            }
            fprintf(stdout, "\n");
            fprintf(stdout, "cov_change (single):\n");
            for (size_t l = 0; l < cov_change.size(); l++) {
                for (size_t k = 0; k < cov_change.at(l).size(); k++) {
                    ss = 0;
                    pp = 0;
                    fprintf(stdout, "|");
                    for (size_t j = 0; j < cov_change[l][k].size(); j++) {
                        ss += cov_change[l][k][j];
                        pp++;
                        fprintf(stdout, "%lu ", cov_change[l][k][j]);
                    }
                    //fprintf(stdout, "%i - %lu|", pp, ss);
                }
                fprintf(stdout, "\n");
            }
            fprintf(stdout, "\n");

            fprintf(stdout, "genome pos:\n");
            for (size_t l = 0; l < genome_pos.size(); l++) {
                for (size_t k = 0; k < genome_pos.at(l).size(); k++) {
                    fprintf(stdout, "|");
                    for (set<unsigned long>::iterator j = genome_pos.at(l).at(k).begin(); j != genome_pos.at(l).at(k).end(); j++)
                        fprintf(stdout, "%lu ", *j);
                }
                fprintf(stdout, "\n");
            }
            fprintf(stdout, "\n");
            
            fprintf(stdout, "cand loss: %f best loss: %f\n\n", candidate_loss, best_loss);
        }
    }

    if (conf->use_mip_objective && ! used_mip) {
        loss = -1.0;
        gain = 0.0;
    } else {
        loss = min(candidate_loss, best_loss);
        gain = max(best_loss - candidate_loss, 0.0);
    }
    return (candidate_loss < best_loss);
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

void get_plifs_from_file() {
    
    FILE* infile = fopen(conf->lossfile.c_str(), "r");    
	if (!infile){
		fprintf(stderr, "Could not open file %s for reading!\n", conf->lossfile.c_str());
		exit(2);
	}
    char* ret;
    char line[10000];
    int i_count;
    float c, sl1, sl2, sr1, sr2;

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
        tmp_plif.push_back((double) sl1);
        tmp_plif.push_back((double) sl2);
        tmp_plif.push_back((double) sr1);
        tmp_plif.push_back((double) sr2);
        plifs.insert(plifs.begin(), pair<double, vector<double> >((double) c, tmp_plif));
    }

    genData->plifs = plifs;
    fclose(infile);
}


double compute_mip_loss(double observed_cov, double predicted_cov, unsigned long segment_len) {

    // in the case of exonic segments, the observed coverage is a count value and the
    // predicted coverage is a mean that needs to be transfered into counts
    // in the case of intronic segments, the coverages are already counts, segment_len is 0 then 
    if (segment_len > 0 && conf->read_len > 0) {
        predicted_cov *= ((double) segment_len / (double) conf->read_len); 
        observed_cov /= (double) conf->read_len; 
    }

    if (observed_cov > 30000)
        observed_cov = 30000;
    if (observed_cov < 0)
        observed_cov = 0;
    
    // get plif iterators to interpolate entries
    map< double, vector<double> >::iterator upper = genData->plifs.lower_bound(observed_cov);
    map< double, vector<double> >::iterator lower = upper;
    if (observed_cov > 0)
        lower--;

    double loss = 0.0;
    double diff = abs(predicted_cov - observed_cov);

    // determine if left or right part of function needs to be queried
    if (predicted_cov > observed_cov) {
        loss = ((lower->second.at(0) + upper->second.at(0)) / 2)*diff*diff + ((lower->second.at(1) + upper->second.at(1)) / 2)*diff; 
    } else {
        loss = ((lower->second.at(2) + upper->second.at(2)) / 2)*diff*diff + ((lower->second.at(3) + upper->second.at(3)) / 2)*diff; 
    }

    return loss;
}


void prepare_mip_objective() {

        if (conf->segmentfile.size() == 0) {
            fprintf(stderr, "For using the mip objective, you need to specify a segments file via -s/--segmentfile. Bailing out!\n");
            exit(-1);
        } else {
            if (conf->verbose) fprintf(stdout, "Parsing segments from %s ...\n", conf->segmentfile.c_str());
            genData->segments.get_from_file(); 
            if (conf->verbose) fprintf(stdout, "... done.\n\n");
        }

       // fprintf(stdout, "chr_size:%i\tchr_num:%i\tcov_map:%i\n", genData->chr_num.size(), genData->chr_num.size(), genData->coverage_map.size());

        if (conf->lossfile.size() == 0) {
            fprintf(stderr, "For using the mip objective, you need to specify a loss function parameter file via -l/--lossfile. Bailing out!\n");
            exit(-1);
        } else {
            if (conf->verbose) fprintf(stdout, "Parsing loss function parameters from %s ...\n", conf->lossfile.c_str());
            get_plifs_from_file();
            if (conf->verbose) fprintf(stdout, "... done.\n\n");
        }
}


void add_zero_segments() {
    
    // handle missing exonic segments
    if (conf->verbose)
        fprintf(stdout, "Adding additional exonic segments for covered but not predicted regions ...\n");

    for (map <pair<unsigned int, unsigned char>, vector<unsigned int> >::iterator it = genData->coverage_map.begin(); it != genData->coverage_map.end(); it++) {
        if (conf->verbose)
            fprintf(stdout, "   ... processing chr %i / %c\n", it->first.first, it->first.second);
        // get boolean vector of covered positions
        vector<bool> coverage(it->second.size(), false);
        for (vector<unsigned int>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
            coverage.at(distance(it->second.begin(), it2)) = (*it2 > 0);
        }
        // set coverage map to false, where we have predicted segments
        if (genData->segments.exons.find(it->first) != genData->segments.exons.end() && genData->segments.exons[it->first].size() > 1) {
            
            set<long> processed_ids;
            for (map<long, long>::iterator it2 = genData->segments.exons[it->first].begin(); it2 != genData->segments.exons[it->first].end(); it2++) {
                if (processed_ids.find(it2->second) != processed_ids.end())
                    continue;
                unsigned long start = genData->segments.exon_ids[it2->second]->start;
                unsigned long end = start + genData->segments.exon_ids[it2->second]->length;

                for (unsigned long i = start; i < end; i++) {
                    coverage.at(i) = false;
                }
                processed_ids.insert(it2->second);
            }
        }

        // infer segments from remaining covered positions
        unsigned int segment_counter = 0;
        vector<bool>::iterator uncov_start = coverage.begin();
        bool in_segment = false;
        for (vector<bool>::iterator c = coverage.begin(); c != coverage.end(); c++) {
            if (*c && in_segment)
                continue;
            else if (*c && ! in_segment) {
                in_segment = true;
                uncov_start = c;
            }
            else if (!*c && ! in_segment)
                continue;
            else {
                long start = distance(coverage.begin(), uncov_start);
                long length = distance(uncov_start, c) + 1;
                Segment* seg = new Segment(start, length, it->first.first, it->first.second, 0.0);
                long segment_id = genData->segments.exon_ids.size() + genData->segments.introns_by_ids.size();
                if (genData->segments.exons.find(it->first) != genData->segments.exons.end()) {
                    genData->segments.exons[it->first].insert(pair<long, long>(start, segment_id));
                    genData->segments.exons[it->first].insert(pair<long, long>(start + length - 1, segment_id));
                } else {
                    map<long, long> tmp_map;
                    tmp_map.insert(pair<long, long>(start, segment_id));
                    tmp_map.insert(pair<long, long>(start + length - 1, segment_id));
                    genData->segments.exons.insert(pair<pair<unsigned char, unsigned char>, map<long, long> >(it->first, tmp_map));
                }
                //fprintf(stdout, "start: %i  end:%i \n", start, start + length - 1);
                genData->segments.exon_ids.insert(pair<long, Segment*>(segment_id, seg));
                segment_counter++;
                in_segment = false;
            }
        }
        if (conf->verbose)
            fprintf(stdout, "   ... added %i segments\n", segment_counter);
    }
    if (conf->verbose)
        fprintf(stdout, "... done.\n\n");

    // handle missing intronic segments
    if (conf->verbose)
        fprintf(stdout, "Adding additional intronic segments for covered but not predicted regions ...\n");

    // iterate over introns found in the alignment file
    for (map <pair<unsigned int, unsigned char>, map< pair<unsigned long, unsigned long>, unsigned int> >::iterator curr_chr = genData->intron_coverage_map.begin(); curr_chr != genData->intron_coverage_map.end(); curr_chr++) {
        if (conf->verbose)
            fprintf(stdout, "   ... processing chr %i / %c\n", curr_chr->first.first, curr_chr->first.second);

        unsigned int segment_counter = 0;
        // check if we have annotated segments for current chromosome
        if (genData->segments.introns.find(curr_chr->first) != genData->segments.introns.end()) {
            // iterate over all alignment-introns of current chr/strand-pair
            for (map< pair<unsigned long, unsigned long>, unsigned int>::iterator aln_int = curr_chr->second.begin(); aln_int != curr_chr->second.end(); aln_int++) {
                // the following gives the same result as equal_range, but is a bit more clear to read
                multimap<unsigned long, unsigned long>::iterator range_pos = genData->segments.introns[curr_chr->first].lower_bound(aln_int->first.first);
                multimap<unsigned long, unsigned long>::iterator range_end = genData->segments.introns[curr_chr->first].upper_bound(aln_int->first.first);
                bool found = false;
                for (; range_pos != range_end; range_pos++) {
                    assert(genData->segments.introns_by_ids[range_pos->second]->start == range_pos->first);
                    // match coordinates of aln_intron and segment matches
                    if (genData->segments.introns_by_ids[range_pos->second]->start + genData->segments.introns_by_ids[range_pos->second]->length - 1 == aln_int->first.second) {
                        found = true;
                        break;
                    }
                }
                // add new intronic segment if we found no segment match
                if (!found) {
                    long intron_id = genData->segments.exon_ids.size() + genData->segments.introns_by_ids.size();
                    Segment* seg = new Segment(aln_int->first.first, aln_int->first.second - aln_int->first.first + 1, curr_chr->first.first, curr_chr->first.second, 0.0);
                  //  fprintf(stdout, "added new intron from %i to %i\n", aln_int->first.first, aln_int->first.second); 
                    genData->segments.introns[curr_chr->first].insert(pair<long, long>(aln_int->first.first, intron_id));
                    genData->segments.introns_by_ids.insert(pair<long, Segment*>(intron_id, seg));
                    segment_counter++;
                }
            }
            if (conf->verbose)
                fprintf(stdout, "   ... added %i segments\n", segment_counter);
        } else {
            // there are no segments for the current chromosome --> add all
            for (map< pair<unsigned long, unsigned long>, unsigned int>::iterator aln_int = curr_chr->second.begin(); aln_int != curr_chr->second.end(); aln_int++) {
                long intron_id = genData->segments.exon_ids.size() + genData->segments.introns_by_ids.size();
                Segment* seg = new Segment(aln_int->first.first, aln_int->first.second - aln_int->first.first + 1, curr_chr->first.first, curr_chr->first.second, 0.0);
                //fprintf(stdout, "added new intron from %i to %i\n", aln_int->first.first, aln_int->first.second); 
                if (genData->segments.introns.find(curr_chr->first) != genData->segments.introns.end()) {
                    genData->segments.introns[curr_chr->first].insert(pair<long, long>(aln_int->first.first, intron_id));
                } else {
                    multimap<unsigned long, unsigned long> tmp_map;
                    tmp_map.insert(pair<unsigned long, unsigned long>(aln_int->first.first, intron_id));
                    genData->segments.introns.insert(pair<pair<unsigned char, unsigned char>, multimap<unsigned long, unsigned long> >(curr_chr->first, tmp_map));
                }
                genData->segments.introns_by_ids.insert(pair<long, Segment*>(intron_id, seg));
                segment_counter++;
            }
        }
    }
    if (conf->verbose)
        fprintf(stdout, "... done.\n\n");
}


bool pair_is_valid(vector<Alignment>::iterator align1, vector<Alignment>::iterator align2) {

    // check for same chromosome, same strand and opposite orientation 
    if ( (align1->chr != align2->chr) || (align1->strand != align2->strand) || (align1->reversed == align2->reversed) )
        return false;
    
    // check for fragment size limit
    int frag_size = max( abs( (int) align1->get_end() - (int) align2->start), abs( (int) align2->get_end() - (int) align1->start));
    if (frag_size >= conf->max_gen_frag_size)
        return false;

    return true;
}


void parse_annotation() {

    if (conf->verbose)
        fprintf(stderr, "\nParsing segment boundaries from annotation file: %s\n", conf->annotation.c_str());

    FILE* infile = fopen(conf->annotation.c_str(), "r");
	if (!infile){
		fprintf(stderr, "Could not open annotation file %s for reading!\n", conf->segmentfile.c_str());
		exit(2);
	}

    char* ret;
    char line[10000];

    while (true) {
        ret = fgets(line, sizeof(line), infile);

        if (!ret)
            break;

        if (line[0] == '#')
            continue;

        char* sl = strtok(line, "\t");
        int idx = 0;
        unsigned int chr = 1;
        char strand = '+';
        long start = 0;
        long stop = 0;
        bool is_exon = false;
        string chr_name;

        while (sl != NULL) {
            if (idx == 0) { // contig
                chr_name = string(sl);
                chr = genData->chr_num[chr_name];
                if (chr == 0) {
                    fprintf(stderr, "\nERROR: The contig names in annotation file seem not to match the ones given in the alignment header. Could not find contig: %s\n", sl);
                    exit(2);
                }
            } else if (idx == 2) { // is exon?
                is_exon = strcmp("exon", sl) == 0;
                if (!is_exon)
                    break;
            } else if (idx == 3) { // start - assumes 1 based inclusive intervals as defined by ensembl GTF
                start = atoi(sl) - 1;
            } else if (idx == 4) { // stop - assumes 1 based incusive intervals as defined by ensembl GTF
                stop = atoi(sl) - 1;
            } else if (idx == 6 && conf->strand_specific) { // strand
                strand = *sl;
            } else if (idx > 4) { // ignore rest
                break;
            }
            sl = strtok(NULL, "\t");
            idx++;
        }

        // only care about exonic segments
        if (is_exon) {
            pair<unsigned int, unsigned char> chr_strand = pair<unsigned int, unsigned char>(chr, strand);
            if (genData->breakpoint_map.find(chr_strand) != genData->breakpoint_map.end()) {
                // annotation exceeds contig in BAM
                if (stop > (long) genData->breakpoint_map[chr_strand].size()) {
                    fprintf(stderr, "WARNING: Segment end %li in contig %s provided in %s is outside the contig length provided in the alignment header (%lu).\n", stop, chr_name.c_str(), conf->annotation.c_str(), genData->breakpoint_map[chr_strand].size());
                    fprintf(stderr, "\tIGNORING SEGMENT\n");
                // probably just a +/- 1 interval convention error
                } else if (stop == (long) genData->breakpoint_map[chr_strand].size()) {
                    genData->breakpoint_map[chr_strand].at(start) = true;
                    genData->breakpoint_map[chr_strand].at(stop - 1) = true;
                // all fine
                } else {
                    genData->breakpoint_map[chr_strand].at(start) = true;
                    genData->breakpoint_map[chr_strand].at(stop) = true;
                }
            }
        }
    }
}

void check_sorted_input() {

    char line[10000];

    FILE* infile = open_bam_pipe_in(conf->infile);

    if (! infile) {
        fprintf(stderr, "Could not open %s for reading!\n", conf->infile.c_str());
        exit(1);
    }
    char* ret = fgets(line, 10000, infile);
    unsigned int counter = 0;

    if (!ret) {
        fprintf(stderr, "Could not read SAM file %s\n", conf->infile.c_str());
        exit(1);
    }

    if (conf->verbose)
        fprintf(stdout, "\nChecking input file %s\n", conf->infile.c_str());

    while (line[0] == '@')
        ret = fgets(line, sizeof(line), infile);

    char* sl = strtok(line, "\t");
    string last_id = string(sl);

    ret = fgets(line, sizeof(line), infile);
    sl = strtok(line, "\t");

    while(ret && (counter < 1000)) {
        string curr_id = string(sl);
        if (strnum_cmp(curr_id.c_str(), last_id.c_str()) < 0) {
            fprintf(stderr, "ERROR: Input file %s seems unsorted:\n\n", conf->infile.c_str());
            fprintf(stderr, "\t ID in record %i: %s\n\t ID in record %i: %s\n\n", counter+1, last_id.c_str(), counter+2, curr_id.c_str());
            fprintf(stderr, "MMR expects input to be sorted by read ID. Please use\n\t samtools sort -n <bamfile>\nto sort your input.\n\n");
            fprintf(stderr, "If your input is sorted by read ID but does not follow the samtools sort convention,\nyou can disable this check using --no-sort-check\n");

            exit(1);
        }
        last_id = curr_id;
        ret = fgets(line, sizeof(line), infile);
        sl = strtok(line, "\t");
        counter++;
    }
    if (conf->verbose) 
        fprintf(stdout, "\tDone - File correctly sorted.\n");
}

