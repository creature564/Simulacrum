

#ifndef sequencing_data_h
#define sequencing_data_h

#include <stdbool.h>
#include <stdint.h>


typedef struct sequencing_data SEQDATA;

extern SEQDATA *newSEQDATA(uint32_t, uint32_t);
extern void freeSEQDATA(SEQDATA *);
extern void displaySEQDATA(SEQDATA *);

extern uint32_t numRegions(SEQDATA *);
extern uint8_t numKmers(SEQDATA *);

extern uint32_t get_region_id(SEQDATA *, const char *);
extern char *get_region_name(SEQDATA *, uint32_t);
extern uint32_t get_region_length(SEQDATA *, uint32_t);
extern float get_region_cov(SEQDATA *, uint32_t);
extern float get_region_gc(SEQDATA *, uint32_t);
extern bool get_region_blast(SEQDATA *, uint32_t);
extern bool get_region_tax(SEQDATA *, uint32_t);
extern char *get_region_histfile(SEQDATA *, uint32_t, uint32_t);

extern void update_region_name(SEQDATA *, uint32_t, const char *);
extern void update_region_length(SEQDATA *, uint32_t, uint32_t);
extern void update_region_cov(SEQDATA *, uint32_t, uint64_t);
extern void update_region_gc(SEQDATA *, uint32_t, uint64_t);
extern void update_region_blast(SEQDATA *, uint32_t, bool);
extern void update_region_tax(SEQDATA *, uint32_t, bool);
extern void update_region_histfile(SEQDATA *, uint32_t, uint32_t, char *);


#endif
