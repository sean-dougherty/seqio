#pragma once

#include "seqio.h"

#include <cstdio>
#include <string>

namespace seqio {
    namespace impl {

        class exception {
        public:
            seqio_err_info const err_info;

            exception(seqio_err_info err_info_) 
                : err_info(err_info_) {
            }
        };

        class sequence_iterator_if {
        public:
            virtual ~sequence_iterator_if() {}

            virtual bool has_next_sequence() = 0;
            virtual seqio_sequence *next_sequence() = 0;
        };

        class metadata_if {
        public:
            virtual ~metadata_if() {}

            virtual uint32_t get_key_count() = 0;
            virtual char const *get_key(uint32_t key_index) = 0;
            virtual char const *get_value(char const *key) = 0;
        };

        class sequence_if {
        public:
            virtual ~sequence_if() {}

            virtual metadata_if *get_metadata() = 0;
            virtual uint32_t read(char *buffer,
                                  uint32_t buffer_length) = 0;
        };

    }
}

#define raise(STATUS, MSG...) {                                         \
            char message[4096];                                         \
            sprintf(message, MSG);                                      \
            seqio_err_info err_info = {SEQIO_ERR_##STATUS, message};    \
            throw seqio::impl::exception(err_info);                     \
        }
