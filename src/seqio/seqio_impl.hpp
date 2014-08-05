#pragma once

#include "seqio.h"

#include <cstdio>

namespace seqio {
    namespace impl {

        class Exception {
        public:
            seqio_err_info const err_info;

            Exception(seqio_err_info err_info_) 
                : err_info(err_info_) {
            }
        };

        class IMetadata {
        public:
            virtual ~IMetadata() {}

            virtual uint32_t getKeyCount() const = 0;
            virtual char const *getKey(uint32_t key_index) const = 0;
            virtual char const *getValue(char const *key) const = 0;
        };

        class ISequence {
        public:
            virtual ~ISequence() {}

            virtual IMetadata const &getMetadata() = 0;
            virtual uint32_t read(char *buffer,
                                  uint32_t buffer_length) = 0;
        };

        class ISequenceIterator {
        public:
            virtual ~ISequenceIterator() {}

            virtual ISequence *nextSequence() = 0;
        };

    }
}

#define raise(STATUS, MSG...) {                                         \
            char message[4096];                                         \
            sprintf(message, MSG);                                      \
            seqio_err_info err_info = {SEQIO_ERR_##STATUS, message};    \
            throw seqio::impl::Exception(err_info);                     \
        }

#define raise_parm(MSG...) {                    \
        raise(INVALID_PARAMETER, MSG);          \
    }

#define raise_io(MSG...) {                      \
    raise(IO, MSG);                             \
    }

#define raise_oom(MSG...) {                      \
        raise(OUT_OF_MEMORY, MSG);               \
    }

#define implement() {                                                   \
        fprintf(stderr, "implement %s @ %s:%d\n", __FUNCTION__, __FILE__, __LINE__); \
        abort();                                                        \
    }
