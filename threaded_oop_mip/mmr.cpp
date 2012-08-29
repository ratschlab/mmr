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
#include "Alignment.h"
#include "GeneralData.h"
#include "OnlineData.h"
#include "BatchData.h"
#include "Utils.h"

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
pthread_mutex_t mutex_counter;

pthread_cond_t fifo_is_full;
pthread_attr_t attr;

queue<OnlineData*> fifo;

void write_output_direct() {

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
    Alignment curr_alignment;

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
        id = curr_alignment.fill(sl, pair_info);
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

void write_output(unordered_map <string, vector<Alignment> > &read_map_left, unordered_map <string, vector<Alignment> > &read_map_right, map<string, unsigned char> &chr_num) {

    FILE* outfile = open_bam_pipe_out(conf->outfile);
    FILE* infile = open_bam_pipe_in(conf->infile);

    char* ret;
    char line[1000];
    char cp_line[1000];

    unsigned char pair_info = 0;

    int output_counter = 0;

    string id;
    Alignment curr_alignment;

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
        curr_alignment.clear();
        id = curr_alignment.fill(sl, pair_info);
        if (id.size() == 0) {
            continue ;
        }
    
        vector<Alignment>::iterator it;
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
            OnlineData* data;
            data = fifo.front();
            fifo.pop();
            pthread_mutex_unlock (&mutex_fifo);
            pthread_mutex_unlock (&mutex_done);
            pthread_cond_signal (&fifo_is_full);
            data->process_data_online(genData);
        } else {
            pthread_mutex_unlock (&mutex_fifo);
            pthread_mutex_unlock (&mutex_done);
        }
    }
    pthread_exit((void*) 0);
}

int main(int argc, char *argv[]) {

    conf = new Config::Config(argc, argv);

    if (conf->verbose)
        conf->print_call(argv[0]);

    if (conf->parse_complete) {
        BatchData* data = new BatchData();
        data->parse_file();

        // pre filter the alignments
        if (conf->pre_filter) {
            data->pre_filter_alignment_maps();
        }

        if (conf->verbose) 
            fprintf(stdout, "\nBuilding up the coverage map ...\n");

        // build coverage map
        unordered_map<string, vector<Alignment> >::iterator r_idx;
        unordered_map<string, vector<Alignment> >::iterator l_idx;
        vector<Alignment>::iterator curr_best;
        unsigned char max_qual;
        
        // handle left reads
        for (l_idx = data->read_map_left.begin(); l_idx != data->read_map_left.end(); l_idx++) {
            max_qual = 0;
            for (vector<Alignment>::iterator v_idx = l_idx->second.begin(); v_idx != l_idx->second.end(); v_idx++) {
                if (v_idx->quality > max_qual) {
                    curr_best = v_idx;
                    max_qual = curr_best->quality;
                }
            }
            if (max_qual == 0) 
                curr_best = l_idx->second.begin();
            curr_best->is_best = true;
            curr_best->update_coverage_map(1);
        }
        // handle right reads
        for (r_idx = data->read_map_right.begin(); r_idx != data->read_map_right.end(); r_idx++) {
            max_qual = 0;
            for (vector<Alignment>::iterator v_idx = r_idx->second.begin(); v_idx != r_idx->second.end(); v_idx++) {
                if (v_idx->quality > max_qual) {
                    curr_best = v_idx;
                    max_qual = curr_best->quality;
                }
            }
            if (max_qual == 0) 
                curr_best = r_idx->second.begin();
            curr_best->is_best = true;
            curr_best->update_coverage_map(1);
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

        // computing active set of alignment
        // an alignment is active if it is currently the used alignment of a read
        data->get_active_read_set();

        if (conf->use_mip_objective)
            prepare_mip_objective();

        if (conf->verbose) { 
            sum_min_loss = data->get_total_min_loss();
            fprintf(stdout, "\nObjective before smoothing:\n\tchanged single (paired):  %i/%i (%i/%i)\n\tobjective: %lf\n", num_changed_single, num_ambiguous_single, num_changed_paired, num_ambiguous_paired, sum_min_loss); 
        }

        // iterate to smooth coverage map
        for (unsigned int iteration = 0; iteration < conf->iterations; iteration++) {

            if (conf->verbose) 
                fprintf(stdout, "\n\t... processing iteration %i of %i ...\n", iteration + 1, conf->iterations);

            num_ambiguous_single = 0;
            num_ambiguous_paired = 0;

            if (conf->use_pair_info) {
                num_changed_paired = data->smooth_coverage_map_paired( num_ambiguous_paired);
            } 
            num_changed_single = data->smooth_coverage_map_single(num_ambiguous_single);

            if (conf->verbose) { 
                sum_min_loss = data->get_total_min_loss();
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

        delete data;
    } else { // process data by on the fly parsing
        
        pthread_t* threads = NULL;
        unsigned int tid;

        for (unsigned int iteration = 0; iteration < conf->iterations; iteration++) {
            
            conf->iteration = iteration;
            OnlineData* data = new OnlineData();
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

                // check, if we need to prepare the mip objective
                if (conf->use_mip_objective)
                    prepare_mip_objective();
            
            } else {
                while (line[0] == '@') {
                    ret = fgets(line, sizeof(line), infile);
                }
                strcpy(last_line, line);
            }

            genData->total_loss = 0.0;
            genData->num_altered = 0;

            string last_id = string("");
            // fill FIFO
            while ((ret = data->parse_file(infile, last_line, genData, counter)) || data->left_reads.size() > 0 || data->right_reads.size() > 0) {
                
                // check if we use the first pass only to fill coverager map and infer zero predicted segments
                if (iteration > 0 || !conf->zero_unpred) {
                    if (conf->num_threads < 2) {
                        data->process_data_online(genData);
                    } else {
                        pthread_mutex_lock (&mutex_fifo);
                        if (fifo.size() >= conf->max_fifo_size) pthread_cond_wait(&fifo_is_full, &mutex_fifo);
                        fifo.push(data);
                        pthread_mutex_unlock (&mutex_fifo);
                    }
                }
                
                if (!ret)
                    break;
                else if (iteration > 0 || !conf->zero_unpred)
                    data = new OnlineData();
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
                delete[] threads;
            }

            if (iteration == 0 && conf->zero_unpred) {
               add_zero_segments(); 
            }

            if (conf->verbose) { 
                fprintf(stdout, "\nsuccessfully parsed %i lines\n", counter - 1);
                fprintf(stdout, "changed %i alignments\n", genData->num_altered);
                fprintf(stdout, "total objective: %f\n", (float) genData->total_loss);
                if (conf->use_mip_objective)
                    fprintf(stdout, "total objective (2): %f\n", (float) genData->segments.get_total_loss());
            }
        }

        write_output_direct();
        /*map <int, vector<unsigned short> >::iterator it = genData->coverage_map.begin();
        for (it; it != genData->coverage_map.end(); it++) {
            fprintf(stdout, "cov vec:\n");
            for (vector<unsigned short>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                fprintf(stdout, "%i ", (*it2));
            }
            fprintf(stdout, "\n");
        }*/
    }

    delete conf;
    delete genData;

    return 0;
}
