#pragma once

#include <stdint.h>

namespace seqio {
    typedef struct {
        uint64_t sequence_offset;
        uint64_t packed_bases_offset;
        uint32_t bases_count;
        uint8_t shift;
    } __attribute__((__packed__)) seqfragment_t;

    typedef struct {
        uint64_t entries_filepos;
        uint32_t entries_count;
    } __attribute__((__packed__)) metadata_t;

    typedef struct {
        uint32_t key;
        uint32_t value;
    } __attribute__((__packed__)) metadata_entry_t;

    typedef struct {
        uint64_t bases_count;
        uint64_t packed_bases_filepos;
        uint64_t packed_bases_length;
        uint64_t seqfragments_filepos;
        uint64_t seqfragments_count;

        metadata_t metadata;
    } __attribute__((__packed__)) sequence_t;

    typedef struct {
        uint64_t filepos;
        uint32_t length;
    } __attribute__((__packed__)) string_storage_t;

    typedef struct {
        uint64_t sequences_filepos;
        uint64_t sequences_count;
        uint64_t max_seqfragments_count;
        uint64_t max_packed_bases_length;
        metadata_t metadata;
        string_storage_t string_storage;
    } __attribute__((__packed__)) header_t;
}
