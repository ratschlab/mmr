#include <string>
#include <cstdio>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <unistd.h>

#include "Alignment.h"
#include "GeneralData.h"
#include "config.h"
#include "Utils.h" 

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
	string command = conf->samtools + string(" view -bS -o ") + out_fname +  " - 2> /dev/null" + " && echo samtools for writing subprocess terminated successfully";
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

double get_variance(vector<vector<unsigned long> > &exon_coverage) {

    double total_var = 0.0;
    for (size_t i = 0; i < exon_coverage.size(); i++) {
        if (exon_coverage.size() < 2) {
            continue;
        } else {
            double sum = 0.0;
            for ( vector<unsigned long>::iterator it = exon_coverage.at(i).begin(); it != exon_coverage.at(i).end(); it++) {
                sum += (double) *it;
            }
            double mean = sum / (double) exon_coverage.at(i).size();
            sum = 0.0;
            for ( vector<unsigned long>::iterator it = exon_coverage.at(i).begin(); it != exon_coverage.at(i).end(); it++) {
                sum += pow((double) *it - mean, 2.0);
            }
            total_var += (sum / ((double) exon_coverage.at(i).size() - 1.0));  
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

bool compare_pair(vector<Alignment>::iterator candidate_left, vector<Alignment>::iterator candidate_right, vector<Alignment>::iterator best_left, vector<Alignment>::iterator best_right, double &loss, bool debug) {

    bool used_mip = false;
    double candidate_loss = 0.0;
    double best_loss = 0.0;

    // check how objective is to be computed
    // we need segments overlapping at least one of the two candidates
    if (conf->use_mip_objective) {
       
        // compute overlap between candidate alignments
        set<unsigned long> already_covered_pos;
        set<unsigned long> not_covered_pos;
        set<unsigned long> empty_set;
        set<unsigned long> overlap;
        set<unsigned long> genome_pos_cl = candidate_left->get_genome_pos();
        set<unsigned long> genome_pos_cr = candidate_right->get_genome_pos();
        set<unsigned long> genome_pos_bl = best_left->get_genome_pos();
        set<unsigned long> genome_pos_br = best_right->get_genome_pos();

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
        vector<vector<unsigned long> > cov_keep;
        vector<vector<unsigned long> > cov_change;

        // init vectors of candidate and best alignments
        vector<pair<vector<Alignment>::iterator, bool> > aligns;
        aligns.push_back( make_pair(candidate_left, false) );
        aligns.push_back( make_pair(candidate_right, false) );
        aligns.push_back( make_pair(best_left, true) );
        aligns.push_back( make_pair(best_right, true) );

        // compute coverage
        compute_coverage_loss(aligns, cov_keep, cov_change);
        
        best_loss = get_variance(cov_keep);
        candidate_loss = get_variance(cov_change);

    }

    if (debug) {
        fprintf(stdout, "best loss: %f\n", best_loss);
        fprintf(stdout, "candidate loss: %f\n", candidate_loss);
        fprintf(stdout, "delta (best - cand): %f\n", best_loss - candidate_loss);
    }
 
    if (conf->use_mip_objective && ! used_mip) {
        loss = -1.0;
    } else {
        loss = min(candidate_loss, best_loss);
    }
    return (candidate_loss < best_loss);
}

bool compare_align_iter_start(const pair<vector<Alignment>::iterator, bool> &left, const pair<vector<Alignment>::iterator, bool> &right) {
    return left.first->start < right.first->start;
}

void compute_coverage_loss(vector<pair<vector<Alignment>::iterator,bool> > aligns, vector<vector<unsigned long> > &cov_keep, vector<vector<unsigned long> > &cov_change) {

    // sort alignments by starting position
    sort(aligns.begin(), aligns.end(), compare_align_iter_start); 

//    unsigned long first_start = (conf->window_size <= aligns.at(0).first->start)?aligns.at(0).first->start - conf->window_size:0ul;

    // determine genomic position set of all alignments (including windows)
    //set<unsigned long> genome_pos = aligns.at(0).first->get_genome_pos(conf->window_size);
    vector<set<unsigned long> > genome_pos;
    for (size_t i = 0; i < aligns.size(); i++) {
        genome_pos.push_back(aligns.at(i).first->get_genome_pos(conf->window_size));
    }

    // fill coverage maps
    for (size_t i = 0; i < aligns.size(); i++) {
        vector<unsigned long> tmp;
        cov_keep.push_back(tmp);
        aligns.at(i).first->fill_coverage_vector(cov_keep.at(i));
    }

    cov_change = cov_keep;

    // iterate over alignments and alter coverage maps
    for (size_t i = 0; i < aligns.size(); i++) {
        aligns.at(i).first->alter_coverage_vector(cov_change, genome_pos, aligns.at(i).second);
    }
}

bool compare_single(vector<Alignment>::iterator candidate, vector<Alignment>::iterator best, double &loss, bool debug) {

    bool used_mip = false;
    double candidate_loss = 0.0;
    double best_loss = 0.0;

    pair<double, double> loss_candidate, loss_best;

    debug=false;//true;

    // check how objective is to be computed
    // we need segments overlapping at least one of the two candidates
    if (conf->use_mip_objective) {

        // compute overlap between candidate and best alignments
        set<unsigned long> genome_pos_candidate = candidate->get_genome_pos();
        set<unsigned long> genome_pos_best = best->get_genome_pos();
        set<unsigned long> empty_set;

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
        vector<vector<unsigned long> > cov_keep;
        vector<vector<unsigned long> > cov_change;

        // init vectors of candidate and best alignments
        vector<pair<vector<Alignment>::iterator, bool> > aligns;
        aligns.push_back( make_pair(candidate, false) );
        aligns.push_back( make_pair(best, true) );

        // compute coverage
        compute_coverage_loss(aligns, cov_keep, cov_change);
        
        best_loss = get_variance(cov_keep);
        candidate_loss = get_variance(cov_change);
        if (debug) {
            fprintf(stdout, "candidate (single):\n");
            candidate->print();
            fprintf(stdout, "best (single):\n");
            best->print();
            fprintf(stdout, "cov_keep (single):\n");
            fprintf(stdout, "cand loss: %f best loss: %f\n\n", candidate_loss, best_loss);
        }
    }

    if (conf->use_mip_objective && ! used_mip) {
        loss = -1.0;
    } else {
        loss = min(candidate_loss, best_loss);
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


