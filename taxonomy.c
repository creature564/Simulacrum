

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libgen.h>
#include <assert.h>
#include <pthread.h>
#include <ctype.h>

#include "taxonomy.h"



struct taxonomy_container {

    uint16_t capacity;
    uint16_t size;
    uint32_t *tax_ids;

};



TAX *newTAX(uint16_t capacity) {

    TAX *node = malloc(sizeof(TAX));
    assert(node != 0);
    node->capacity = capacity;
    node->size = 0;
    node->tax_ids = calloc(capacity, sizeof(uint32_t));
    assert(node->tax_ids != 0);
    return node;

}


bool checkTAX(TAX *node, uint32_t tax_id) {

    for(uint16_t i = 0; i < node->size; i++) {
        if(node->tax_ids[i] == tax_id) return true;
    }

    return false;

}


void insertTAX(TAX *node, uint32_t tax_id) {

    assert(node != 0);
    if(!checkTAX(node, tax_id)) {
        if(node->size == node->capacity) {
            uint32_t *tmp = realloc(node->tax_ids, (node->capacity * 2) * sizeof(uint32_t));
            assert(tmp != 0);
            node->tax_ids = tmp;
            node->capacity *= 2;
        }
        node->tax_ids[node->size] = tax_id;
        node->size++;
    }

}



void freeTAX(TAX *node) {

    free(node->tax_ids);
    free(node);

}


void displayTAX(TAX *node){

    printf("blast hits:");
    for(uint16_t i = 0; i < node->size; i++) printf(" %d", node->tax_ids[i]);
    printf("\n");
}



uint16_t numHits(TAX *node) {

    return node->size;

}


uint32_t getTaxID(TAX *node, uint32_t i) {

    return node->tax_ids[i];

}




TAX **parse_blast(char *filepath, SEQDATA *data) {

    TAX **blast_table = calloc(numRegions(data), sizeof(TAX*));
    assert(blast_table != 0);

    FILE *fp = fopen(filepath, "r");
    assert(fp != 0);

    size_t buff_size = 256;
    char *line_buff = malloc(buff_size);
    ssize_t read;
    char *prev_region_name = malloc(buff_size);
    uint32_t prev_seq_id = 0, prev_tax_id = 0;

    while ((read = getline(&line_buff, &buff_size, fp)) != -1) {

        char *tok = strchr(line_buff, '\t');
        uint32_t tax_id = (uint32_t)strtol(tok, 0, 10);
        size_t region_name_len = strlen(line_buff) - strlen(tok);
        char region_name[region_name_len];
        strncpy(region_name, line_buff, region_name_len);
        region_name[region_name_len] = '\0';

        if(strcmp(region_name, prev_region_name) == 0) {
            if(tax_id == prev_tax_id) continue;
            else insertTAX(blast_table[prev_seq_id], tax_id);
        }
        else {
            uint32_t seq_id = get_region_id(data, region_name);
            update_region_blast(data, seq_id, true);
            blast_table[seq_id] = newTAX(2);
            insertTAX(blast_table[seq_id], tax_id);
            strcpy(prev_region_name, region_name);
            prev_seq_id = seq_id;
        }

        prev_tax_id = tax_id;

    }

    fclose(fp);
    free(prev_region_name);
    free(line_buff);

    return blast_table;

}



bool check_delnodes(const char *filepath, const uint32_t check_id) {

    FILE *fp = fopen(filepath, "r");
    assert(fp != 0);

    size_t buff_size = 32;
    char *line_buff = malloc(buff_size);
    assert(line_buff != 0);
    ssize_t read;

    while ((read = getline(&line_buff, &buff_size, fp)) != -1) {
        if(check_id == (uint32_t)strtol(line_buff, 0, 10)) return true;
    }

    fclose(fp);
    free(line_buff);

    return false;

}



void check_merged(const char *filepath, uint32_t *id_addr) {

    FILE *fp = fopen(filepath, "r");
    assert(fp != 0);

    size_t buff_size = 64;
    char *line_buff = malloc(buff_size);
    assert(line_buff != 0);
    ssize_t read;

    while ((read = getline(&line_buff, &buff_size, fp)) != -1) {
        char *line_r;
        int32_t curr_id = (uint32_t)strtol(strtok_r(line_buff, "|", &line_r), 0, 10);
        uint32_t new_id = (uint32_t)strtol(strtok_r(NULL, "|", &line_r), 0, 10);
        if(*id_addr == curr_id) {
            memcpy(id_addr, &new_id, sizeof(uint32_t));
            break;
        }
    }

    fclose(fp);
    free(line_buff);

}



void check_nodes(const char *filepath, uint32_t *id_addr, const uint32_t check_id, const char *check_rank) {

    FILE *fp = fopen(filepath, "r");
    assert(fp != 0);

    size_t buff_size = 256;
    char *line_buff = malloc(buff_size);
    assert(line_buff != 0);
    ssize_t read;

    while ((read = getline(&line_buff, &buff_size, fp)) != -1) {
        char *line_r;
        uint32_t curr_id = (uint32_t)strtol(strtok_r(line_buff, "|", &line_r), 0, 10);
        if(check_id == curr_id) {
            uint32_t parent_id = (uint32_t)strtol(strtok_r(NULL, "|", &line_r), 0, 10);
            char *rank_tok = strtok_r(NULL, "|", &line_r);
            if(strstr(rank_tok, check_rank)) {
                memcpy(id_addr, &curr_id, sizeof(uint32_t));
            }
            else {
                if(curr_id == parent_id) break; // root
                else check_nodes(filepath, id_addr, parent_id, check_rank);
            }
            break;
        }
    }

    fclose(fp);
    free(line_buff);

}



bool check_names(const char *filepath, const uint32_t check_id, const char *check_class) {

    FILE *fp = fopen(filepath, "r");
    assert(fp != 0);

    size_t buff_size = 256;
    char *line_buff = malloc(buff_size);
    assert(line_buff != 0);
    ssize_t read;

    while ((read = getline(&line_buff, &buff_size, fp)) != -1) {
        char *line_r;
        uint32_t curr_id = (uint32_t)strtol(strtok_r(line_buff, "|", &line_r), 0, 10);
        char *class = strlower(strtok_r(NULL, "|", &line_r));
        if(check_id == curr_id) {
            if(strstr(line_r, "scientific name")) {
              if(strstr(class, check_class)) return true;
              else break;
            }
            else continue;
        }
    }

    fclose(fp);
    free(line_buff);

    return false;

}


char *strlower(char *str) {

    uint8_t len = strlen(str);
    for(uint8_t i = 0; i < len; i++){
        str[i] = tolower(str[i]);
    }
    return str;
}
