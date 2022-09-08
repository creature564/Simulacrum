


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#include "sequence.h"
#include "hashmap.h"



int compareKNODE(void *c, void *k) {
    KNODE *check = (KNODE*)c;
    KNODE *knode = (KNODE*)k;
    if(check->kmer == knode->kmer) {
        knode->count++;
        return 0;
    }

    else return 1;

}



uint64_t hashKNODE(const void *item, uint64_t seed0, uint64_t seed1) {
    //KNODE *knode = (KNODE *)item;
    return hashmap_murmur(item, sizeof(uint64_t), seed0, seed1);
}



bool iterKNODE(const void *k) {
    KNODE *knode = (KNODE *)k;
    printf("%lu\t%d\n", knode->kmer, knode->count);
    return true;
}



struct sequence_encoding {

    uint32_t capacity;
    uint32_t size;
    uint64_t *code;

};



SEQCODE *newSEQCODE(uint32_t seq_len) {

  SEQCODE *seqcode = malloc(sizeof(SEQCODE));
  assert(seqcode != 0);
  seqcode->capacity = seq_len / 32;
  if(seq_len % 32 != 0) seqcode->capacity++;
  seqcode->size = 0;
  seqcode->code = calloc(seqcode->capacity, sizeof(uint64_t));
  assert(seqcode->code != 0);
  return seqcode;

}



void insertSEQCODE(SEQCODE *seqcode, uint64_t code) {

    assert(seqcode->size < seqcode->capacity);
    seqcode->code[seqcode->size] = code;
    seqcode->size++;

}


uint64_t getSEQCODE(SEQCODE *seqcode, uint32_t i) {

    return seqcode->code[i];

}



void freeSEQCODE(SEQCODE *seqcode) {

    free(seqcode->code);
    free(seqcode);

}



uint64_t encode_count_gc(const char *seq, uint32_t len, uint64_t *gc_count) {

    uint64_t code = 0;

    for(uint32_t i = 0; i < len; i++) {
        uint32_t j = i % 32;
        switch(seq[i]) {
      //      case 'A':
      //      case 'a':
      //          code |= ((uint64_t)0x0 << (j * 2));
      //          break;
              case 'C':
              case 'c':
                  code |= ((uint64_t)0x1 << (j * 2));
                  (*gc_count)++;
                  break;
              case 'G':
              case 'g':
                  code |= ((uint64_t)0x2 << (j * 2));
                  (*gc_count)++;
                  break;
              case 'T':
              case 't':
                  code |= ((uint64_t)0x3 << (j * 2));
                  break;
          }
    }

    //printf("%lu\n", code);
    //printBits(code);
    //printf("\n");
    return code;

}


// zero-based indexing
uint64_t get_kmer(SEQCODE *seqcode, uint8_t k, uint32_t pos) {

    uint64_t kmer, buffer, overflow;
    uint32_t div, mod;

    kmer = 0;
    div = pos / 32;
    mod = pos % 32;
    //printf("mod: %lu\n", mod);
    //printf("k: %lu\n", k);
    buffer = seqcode->code[div];
    //printf("buffer:\n");
    //printBits(buffer);

    if(32 - mod < k) {
        //printf("32 - mod < k\n\n");
        buffer >>= (mod * 2);
        //printf("buffer >>= mod * 2:\n");
        //printBits(buffer);
        kmer |= buffer;
        //printf("kmer |= buffer:\n");
        //printBits(kmer);
        overflow = seqcode->code[div + 1];
        //printf("overflow:\n");
        //printBits(overflow);
        overflow <<= ((32 - mod) * 2);
        //printf("overflow <<= (32 - mod) * 2:\n");
        //printBits(overflow);
        kmer |= (overflow << ((32 - k) * 2) >> ((32 - k) * 2));
        //printf("kmer |= (overflow << ((32 - k) * 2) >> ((32 - k) * 2)):\n");
        //printBits(kmer);
    }
    else if(32 - mod == k) {
        //printf("32 - mod == k\n\n");
        buffer >>= (mod * 2);
        //printf("buffer >>= mod * 2:\n");
        //printBits(buffer);
        kmer |= buffer;
        //printf("kmer |= buffer:\n");
        //printBits(kmer);
    }
    else {
        //printf("32 - mod > k\n\n");
        kmer |= (buffer << ((32 - mod - k) * 2) >> ((32 - k) * 2));
        //printf("kmer |= (buffer << ((32 - mod - k) * 2) >> ((32 - k) * 2)):\n");
        //printBits(kmer);
    }

    return kmer;
}


// https://stackoverflow.com/questions/52845040/printing-a-long-in-binary-64-bit-representation
void printBits(uint64_t n){
    uint64_t i = 1UL<<(sizeof(n)*CHAR_BIT-1);
    while(i>0){
         if(n&i)
              printf("1");
         else
              printf("0");
         i >>= 1;
    }
    printf("\n\n");
}
