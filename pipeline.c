
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <assert.h>

#include "htslib/thread_pool.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "pipeline.h"
#include "taxonomy.h"
#include "sequence.h"
#include "hashmap.h"



static void *pipe_in2tax(void *arg);
static void *process_taxonomy(void *arg);
static void *pipe_tax2seq(void*arg);
static void *process_sequence(void *arg);



typedef struct {
    hts_tpool *p;
    hts_tpool_process *q1;
    hts_tpool_process *q2;
    PARAM *param;
    SEQDATA *seqdata;
    TAX **blast_table;
    TAX *check_table;
    pthread_mutex_t *lock;
} pipe_opt;



typedef struct {
    pipe_opt *o;
    uint32_t id;
} pipe_job;



static void *pipe_in2tax(void *arg) {
    pipe_opt *o = (pipe_opt *)arg;
    for (uint32_t i = 0; i < numRegions(o->seqdata); i++) {
        pipe_job *job = malloc(sizeof(pipe_job));
        assert(job != 0);
        job->o = o;
        job->id = i;
        if (hts_tpool_dispatch(o->p, o->q1, process_taxonomy, job) != 0) {
            free(job);
            pthread_exit((void *)1);
        }
    }
    pthread_exit(0);
}



static void *process_taxonomy(void *arg) {

    pipe_job *job = (pipe_job *)arg;

    if(job->o->blast_table[job->id] != 0) {

      for(uint16_t i = 0; i < numHits(job->o->blast_table[job->id]); i++) {

          uint32_t curr_TaxID = getTaxID(job->o->blast_table[job->id], i);

          pthread_mutex_lock(job->o->lock);
          if(checkTAX(job->o->check_table, curr_TaxID)) {
              pthread_mutex_unlock(job->o->lock);
              update_region_tax(job->o->seqdata, job->id, true);
              break;
          }
          pthread_mutex_unlock(job->o->lock);

          uint32_t check = curr_TaxID;

          if(check_delnodes(job->o->param->delnodes, check)) continue;
          check_merged(job->o->param->merged, &check);
          check_nodes(job->o->param->nodes, &check, check, job->o->param->rank);
          if(check_names(job->o->param->names, check, job->o->param->classification)) {
              pthread_mutex_lock(job->o->lock);
              insertTAX(job->o->check_table, curr_TaxID);
              pthread_mutex_unlock(job->o->lock);
              update_region_tax(job->o->seqdata, job->id, true);
          }
      }

      freeTAX(job->o->blast_table[job->id]);

    }

    return (void *)job;

}



static void *pipe_tax2seq(void *arg) {
    pipe_opt *o = (pipe_opt *)arg;
    hts_tpool_result *result;
    while ((result = hts_tpool_next_result_wait(o->q1))) {
        pipe_job *job = (pipe_job *)hts_tpool_result_data(result);
        hts_tpool_delete_result(result, 0);
        if (hts_tpool_dispatch(job->o->p, job->o->q2, process_sequence, job) != 0)
            pthread_exit((void *)1);
        if (hts_tpool_process_empty(o->q1)) {
            free(o->blast_table);
            freeTAX(o->check_table);
            break;
        }
    }
    pthread_exit(0);
}



static void *process_sequence(void *arg) {

    pipe_job *job = (pipe_job *)arg;

    faidx_t *fa_idx = fai_load(job->o->param->assembly_path);
    assert(fa_idx != 0);
    samFile *sam = sam_open(job->o->param->alignment_path, "r");
    assert(sam != 0);
    hts_idx_t *hts_idx = sam_index_load(sam, job->o->param->alignment_path);
    assert(hts_idx != 0);
    bam_hdr_t *bam_hdr = sam_hdr_read(sam);
    assert(bam_hdr != 0);

    uint32_t seqlen = bam_hdr->target_len[job->id];

    update_region_length(job->o->seqdata, job->id, seqlen);



    // gc content analysis & sequence encoding

    SEQCODE *seqcode = newSEQCODE(seqlen);

    uint64_t gc_count = 0;
    char *seq_buff;
    uint32_t i = 0, buff_size;
    while(i < seqlen - 32) {
        seq_buff = faidx_fetch_seq(fa_idx, get_region_name(job->o->seqdata, job->id), i, i+31, &buff_size);
        insertSEQCODE(seqcode, encode_count_gc(seq_buff, buff_size, &gc_count));
        free(seq_buff);
        i += 32;
    }
    if(seqlen - i > 0) {
        seq_buff = faidx_fetch_seq(fa_idx, get_region_name(job->o->seqdata, job->id), i, seqlen, &buff_size);
        insertSEQCODE(seqcode, encode_count_gc(seq_buff, buff_size, &gc_count));
        free(seq_buff);
    }

    update_region_gc(job->o->seqdata, job->id, gc_count);

    fai_destroy(fa_idx);



    // kmer analysis

    hashmap *hash = hashmap_new(sizeof(KNODE), seqlen, getSEQCODE(seqcode, 0), getSEQCODE(seqcode, 1), hashKNODE, compareKNODE, iterKNODE);

    for(uint32_t k = 0; k < job->o->param->n_kmers; k++) {

        uint32_t kmer = job->o->param->kmer_list[k];

        for(uint32_t i = 0; i < seqlen - kmer; i++) {
            hashmap_set(hash, &(KNODE){.kmer = get_kmer(seqcode, kmer, i), .count = 1});
        }

        char *histFile = exportKHASH(hash, get_region_name(job->o->seqdata, job->id), kmer, job->o->param->max_kmer_freq);
        update_region_histfile(job->o->seqdata, job->id, k, histFile);
        //free(histFile);
        hashmap_clear(hash);

    }

    hashmap_free(hash);
    freeSEQCODE(seqcode);



    // coverage analysis

    hts_itr_t *itr = sam_itr_querys(hts_idx, bam_hdr, bam_hdr->target_name[job->id]);
    assert(itr != 0);
    bam1_t *itr_data = bam_init1();
    assert(itr_data != 0);

    uint64_t bp_count = 0;
    while (sam_itr_next(sam, itr, itr_data) >= 0) {
        bam1_core_t *c = &(itr_data->core);
        if((c->pos + c->l_qseq) < seqlen) bp_count += c->l_qseq;
        else bp_count += (seqlen - c->pos);
    }

    update_region_cov(job->o->seqdata, job->id, bp_count);

    bam_destroy1(itr_data);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bam_hdr);
    hts_idx_destroy(hts_idx);
    sam_close(sam);
    free(job);

}



SEQDATA *run_pipeline(PARAM *param) {

    faidx_t *fa_idx = fai_load(param->assembly_path);
    assert(fa_idx != 0);
    uint32_t n_seq = faidx_nseq(fa_idx);
    SEQDATA *seqdata = newSEQDATA(param->n_kmers, n_seq);
    for(uint32_t i = 0; i < n_seq; i++) update_region_name(seqdata, i, faidx_iseq(fa_idx, i));
    fai_destroy(fa_idx);

    TAX **blast_table = parse_blast(param->blast_path, seqdata);
    TAX *taxdump_table = newTAX(10);
    pthread_mutex_t lock;
    pthread_mutex_init(&lock, NULL);

    hts_tpool *p = hts_tpool_init(param->threads_count);
    hts_tpool_process *q1 = hts_tpool_process_init(p, n_seq*2, 0);  // taxonomy queue
    hts_tpool_process *q2 = hts_tpool_process_init(p, n_seq, 1);    // sequence queue
    pipe_opt o = {p, q1, q2, param, seqdata, blast_table, taxdump_table, &lock};
    pthread_t tidIto1, tid1to2;
    void *retv;
    int ret;

  // Launch our data source and sink threads.
    pthread_create(&tidIto1, NULL, pipe_in2tax, &o);
    pthread_create(&tid1to2, NULL, pipe_tax2seq, &o);

  // Wait for tasks to finish.
    ret = 0;
    pthread_join(tidIto1, &retv); ret |= (retv != NULL);
    pthread_join(tid1to2, &retv); ret |= (retv != NULL);
    printf("Return value %d\n", ret);

    hts_tpool_process_destroy(q1);
    hts_tpool_process_flush(q2);
    hts_tpool_process_destroy(q2);
    hts_tpool_destroy(p);

    pthread_mutex_destroy(&lock);

    return seqdata;

}
