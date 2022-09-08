
#ifndef sequence_h
#define sequence_h

#include <stdint.h>
#include <stdbool.h>

typedef struct kmer_node {

    uint64_t kmer;
    uint16_t count;

} KNODE;

extern KNODE *newKNODE(uint64_t);
extern int compareKNODE(void *, void *);
extern uint64_t hashKNODE(const void *, uint64_t, uint64_t);
extern bool iterKNODE(const void *);

typedef struct sequence_encoding SEQCODE;

extern SEQCODE *newSEQCODE(uint32_t);
extern void insertSEQCODE(SEQCODE *, uint64_t);
extern uint64_t getSEQCODE(SEQCODE *, uint32_t);
extern void freeSEQCODE(SEQCODE *);
extern uint64_t encode_count_gc(const char *, uint32_t, uint64_t *);
extern uint64_t get_kmer(SEQCODE *, uint8_t, uint32_t);

extern void printBits(uint64_t);

#endif
