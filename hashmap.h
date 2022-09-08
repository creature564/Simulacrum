// Copyright 2020 Joshua J Baker. All rights reserved.
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.

#ifndef HASHMAP_H
#define HASHMAP_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct hashmap hashmap;

struct hashmap *hashmap_new(size_t elsize, size_t cap,
                            uint64_t seed0, uint64_t seed1,
                            uint64_t (*hash)(const void *,
                                             uint64_t, uint64_t),
                            int (*compare)(void *, void *), bool (*iter)(const void *));
void hashmap_free(struct hashmap *map);
void hashmap_clear(struct hashmap *map);
size_t hashmap_count(struct hashmap *map);
bool hashmap_oom(struct hashmap *map);
int hashmap_set(struct hashmap *map, void *);
bool hashmap_scan(struct hashmap *map);
void hashmap_set_allocator(void *(*malloc)(size_t), void (*free)(void*));
uint64_t hashmap_sip(const void*, size_t len, uint64_t seed0, uint64_t seed1);
uint64_t hashmap_murmur(const void*, size_t len, uint64_t seed0, uint64_t seed1);

extern char *exportKHASH(hashmap *, char *, uint32_t, uint32_t);

#endif
