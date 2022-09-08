

#ifndef param_h
#define param_h

#include <stdint.h>

typedef struct {

    uint32_t threads_count;

    char *assembly_path;
    char *alignment_path;
    char *blast_path;

    char *delnodes;
    char *merged;
    char *nodes;
    char *names;

    char *rank;
    char *classification;

    uint32_t max_kmer_freq;
    uint32_t n_kmers;
    uint32_t *kmer_list;

} PARAM;

extern void initPARAM(PARAM *, int, char **);
extern void freePARAM(PARAM *);

#endif
