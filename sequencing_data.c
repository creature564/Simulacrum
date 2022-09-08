

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "sequencing_data.h"



struct sequencing_data {

    uint32_t n_regions;                   // number of regions
    char **names;                         // region names
    uint32_t *lengths;                    // region lengths
    float *cov;                           // alignment coverages per region
    float *gc;                            // gc content % per region
    bool *blast_hits;                     // does region have blast hit?
    bool *target_tax;                     // does blast hit match target taxonomy at desired classification level?
    uint32_t n_kmers;                      // number of kmers to be used in analysis
    char ***hist_files;                   // kmer distribution filenames per region

};



SEQDATA *newSEQDATA(uint32_t n_kmers, uint32_t n_regions) {

    SEQDATA *data = malloc(sizeof(SEQDATA));
    assert(data != 0);

    data->n_kmers = n_kmers;
    data->n_regions = n_regions;

    data->names = calloc(n_regions, sizeof(char*));
    assert(data->names != 0);

    data->lengths = calloc(n_regions, sizeof(uint32_t));
    assert(data->lengths != 0);

    data->cov = calloc(n_regions, sizeof(float));
    assert(data->cov != 0);

    data->gc = calloc(n_regions, sizeof(float));
    assert(data->gc != 0);

    data->blast_hits = calloc(n_regions, sizeof(bool));
    assert(data->blast_hits != 0);

    data->target_tax = calloc(n_regions, sizeof(bool));
    assert(data->target_tax != 0);

    data->hist_files = calloc(n_regions, sizeof(char**));
    assert(data->hist_files != 0);

    for(uint32_t i = 0; i < n_regions; i++) {

        data->blast_hits[i] = false;
        data->target_tax[i] = false;

        data->hist_files[i] = calloc(n_kmers, sizeof(char*));
        assert(data->hist_files[i] != 0);

    }

    return data;

}



void freeSEQDATA(SEQDATA *data) {

    if (data == 0) return;

    if (data->names) {

        for (uint32_t i = 0; i < data->n_regions; i++) {

          free(data->names[i]);
          for(uint32_t k = 0; k < data->n_kmers; k++) free(data->hist_files[i][k]);
          free(data->hist_files[i]);

        }

        free(data->names);
        free(data->hist_files);
        free(data->lengths);
        free(data->cov);
        free(data->gc);
        free(data->blast_hits);
        free(data->target_tax);

    }

    free(data);

}



void displaySEQDATA(SEQDATA *data) {

    for(uint32_t id = 0; id < data->n_regions; id++) {

        printf("Region: %s\n", data->names[id]);
        printf("Length: %d\n", data->lengths[id]);
        printf("GC Content: %f\n", data->gc[id]);
        printf("Read Coverage: %f\n", data->cov[id]);
        printf("BLAST Hit: %s\n", data->blast_hits[id] ? "true" : "false");
        printf("Target taxon: %s\n", data->target_tax[id] ? "true" : "false");
        for(uint32_t k = 0; k < data->n_kmers; k++) printf("%s\n", data->hist_files[id][k]);
        printf("\n");

    }

}



uint32_t numRegions(SEQDATA *data) {

    return data->n_regions;

}



uint8_t numKmers(SEQDATA *data) {

    return data->n_kmers;

}



char *get_region_name(SEQDATA *data, uint32_t id) {

    return data->names[id];

}



uint32_t get_region_length(SEQDATA *data, uint32_t id) {

    return data->lengths[id];

}



float get_region_cov(SEQDATA *data, uint32_t id) {

    return data->cov[id];

}



float get_region_gc(SEQDATA *data, uint32_t id) {

    return data->gc[id];

}



bool get_region_blast(SEQDATA *data, uint32_t id) {

    return data->blast_hits[id];

}



bool get_region_tax(SEQDATA *data, uint32_t id) {

    return data->target_tax[id];

}



char *get_region_histfile(SEQDATA *data, uint32_t r_id, uint32_t k_id) {

    return data->hist_files[r_id][k_id];

}



uint32_t get_region_id(SEQDATA *data, const char *name) {

    for(uint32_t i = 0; i < data->n_regions; i++) {

        if((!data->names[i]) || (strcmp(data->names[i], name) == 0)) return i;

    }

    return -1;

}



void update_region_name(SEQDATA *data, uint32_t id, const char *name) {

    if(!data->names[id]) {

      data->names[id] = malloc(strlen(name+1));
      strcpy(data->names[id], name);

    }

    else assert(strcmp(data->names[id], name) == 0);

}



void update_region_length(SEQDATA *data, uint32_t id, uint32_t length) {

    if(!data->lengths[id]) data->lengths[id] = length;
    else assert(data->lengths[id] == length);

}



void update_region_cov(SEQDATA *data, uint32_t id, uint64_t bp_count) {

    assert(data->lengths[id] != 0);
    float cov = (double)bp_count / (float)data->lengths[id];
    if(!data->cov[id]) data->cov[id] = cov;
    else assert(data->cov[id] == cov);

}



void update_region_gc(SEQDATA *data, uint32_t id, uint64_t gc_count) {

    assert(data->lengths[id] != 0);
    float gc = ((double)gc_count / (float)data->lengths[id]) * 100;
    if(!data->gc[id]) data->gc[id] = gc;
    else assert(data->gc[id] == gc);

}



void update_region_blast(SEQDATA *data, uint32_t id, bool has_blast) {

    if(!data->blast_hits[id]) data->blast_hits[id] = has_blast;
    else assert(data->blast_hits[id] == has_blast);

}



void update_region_tax(SEQDATA *data, uint32_t id, bool is_target) {

    if(!data->target_tax[id]) data->target_tax[id] = is_target;
    else assert(data->target_tax[id] == is_target);

}



void update_region_histfile(SEQDATA *data, uint32_t r_id, uint32_t k_id, char *file) {

    if(!data->hist_files[r_id][k_id]) {

      data->hist_files[r_id][k_id] = file; //malloc(strlen(file+1));
      //strcpy(data->hist_files[r_id][k_id], file);

    }

    else assert(strcmp(data->hist_files[r_id][k_id], file) == 0);

}
