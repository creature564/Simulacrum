cmake_minimum_required(VERSION 3.0.2)
project(SIDR2.0 C)

set(CMAKE_C_STANDARD 99)

message(STATUS "Build: Release")

set(SRC
        htslib/libhts.a
        src/assert.h
        src/logging.h src/logging.c
        src/mem_allocator.h src/mem_allocator.c
        src/queue.h src/queue.c
        src/buffer_queue.h src/buffer_queue.c
        src/file_reader.h src/file_reader.c
        src/file_writer.h src/file_writer.c
        src/hash_map.h src/hash_map.c
        src/kmer_processor.h src/kmer_processor.c
        src/kc_param.h src/kc_param.c
        src/kc_utils.h src/kc_utils.c
        src/kmer_counter.h src/kmer_counter.c
        src/header.h src/header.c
        src/histo.h src/histo.c
        src/pthread_barrier.h src/pthread_barrier.c
        src/chtkc.h src/chtkc.c
        src/param.h src/param.c
        src/pylink.h src/pylink.c
        src/sequence.h src/sequence.c
        src/sequencing_data.h src/sequencing_data.c
        src/taxonomy.h src/taxonomy.c
        src/utils.h src/utils.c
        src/main.c
        )

add_executable(SIDR2.0 ${SRC})
target_link_libraries(SIDR2.0 htslib pthread z m bz2 lzma curl crypto)
target_include_directories(SIDR2.0 htslib)

install(TARGETS SIDR2.0 DESTINATION .)
