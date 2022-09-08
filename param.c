

#include <argp.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "param.h"


static char *getFilepathFromDir(const char *dir, const char *filename) {

    int dlen = strlen(dir);
    int flen = strlen(filename);
    char *filepath = malloc(dlen + 1 + flen + 1);
    assert(filepath != 0);
    if(dir[dlen-1] == '/') strcat(strcpy(filepath, dir), filename);
    else strcat(strcat(strcpy(filepath, dir), "/"), filename);

    return filepath;

}

/*
static void parse_kmers(struct argp_state *state, const char *arg, const char *info) {

    PARAM *param = state->input;

}
*/


static error_t parse_opt(int key, char *arg, struct argp_state *state) {

    PARAM *param = state->input;

    switch (key) {

        case 'a':
            param->assembly_path = arg;
            break;
        case 'b':
            param->blast_path = arg;
            break;
        case 't':
            param->delnodes = getFilepathFromDir(arg, "delnodes.dmp");
            param->merged = getFilepathFromDir(arg, "merged.dmp");
            param->nodes = getFilepathFromDir(arg, "nodes.dmp");
            param->names = getFilepathFromDir(arg, "names.dmp");
            break;
        case 'r':
            param->rank = arg;
            break;
        case 'c':
            param->classification = arg;
            break;
        //case 'k':
        //    parse_kmer(state, arg, "Kmer list");
        //    break;
        case 'n':
            param->threads_count = (uint32_t)strtol(arg, 0, 0);
            break;
        case ARGP_KEY_ARG:
            param->alignment_path = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 1)
                argp_usage(state);
            if (!param->assembly_path)
                argp_error(state, "Assembly file (.fasta/.fa/.fna) must be provided");
            if (!param->alignment_path)
                argp_error(state, "Alignment file (.sam/.bam/.cram) must be provided");
            if (!param->blast_path)
                argp_error(state, "Blast query output file must be provided");
            if (!param->delnodes || !param->merged || !param->nodes || !param->names)
                argp_error(state, "NCBI Taxdump directory must be provided");
            if (!param->classification)
                argp_error(state, "Target classification must be provided");
            break;
        default:
            return ARGP_ERR_UNKNOWN;

    }

    return 0;

}



void initPARAM(PARAM *param, int argc, char **argv) {

    printf("Begin param initialization\n");

    // default values

    param->rank = "phylum";
    param->max_kmer_freq = 256;
    param->n_kmers = 11;
    uint32_t kmer_list[11] = {15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
    param->kmer_list = malloc(sizeof(kmer_list));
    memcpy(param->kmer_list, kmer_list, sizeof(kmer_list));
    param->threads_count = 4;

    struct argp_option options[] = {

            {"fa", 'a', "FILE", 0, "Assembly file (fasta/fa/fna)", 0},
            {"blast", 'b', "FILE", 0, "NCBI Blast+ query file", 0},
            {"tax", 't', "DIRECTORY", 0, "NCBI Taxdump directory", 0},
            {"rank", 'r', "STRING", 0, "Target taxonomic rank, default = 'phylum'", 0},
            {"class", 'c', "STRING", 0, "Target classification to identify", 0},

            {"kmers", 'k', "[a,c-m,p ...]", 0, "Kmers to include in analysis, default = [15-25]", 1},
            {"threads", 'n', "N", 0, "Number of processing threads, default = 4", 1},

            {0}

    };

    struct argp argp = {options, parse_opt, "ARGS...", "Sequence Identification"};
    argp_parse(&argp, argc, argv, 0, 0, param);

    printf("Alignment: %s\n", param->alignment_path);
    printf("Assembly: %s\n", param->assembly_path);
    printf("Blast: %s\n", param->blast_path);
    printf("Delnodes: %s\n", param->delnodes);
    printf("Merged: %s\n", param->merged);
    printf("Nodes: %s\n", param->nodes);
    printf("Names: %s\n", param->names);
    printf("Rank: %s\n", param->rank);
    printf("Classification: %s\n", param->classification);

}

void freePARAM(PARAM *param) {

    free(param->delnodes);
    free(param->merged);
    free(param->nodes);
    free(param->names);
    free(param->kmer_list);

}
