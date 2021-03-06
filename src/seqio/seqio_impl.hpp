#pragma once

#include "seqio.h"

#include <cstdio>
#include <cstdlib>

//#define SEQIO_ABORT_ON_EXCEPTION

namespace seqio {
    namespace impl {

        class Exception {
        public:
            seqio_err_info const err_info;

            Exception(seqio_err_info err_info_) 
                : err_info(err_info_) {

#ifdef SEQIO_ABORT_ON_EXCEPTION
                abort();
#endif
            }
        };

        class IConstDictionary {
        public:
            virtual ~IConstDictionary() {}

            virtual bool hasKey(char const *key) const = 0;
            virtual uint32_t getKeyCount() const = 0;
            virtual char const *getKey(uint32_t key_index) const = 0;
            virtual char const *getValue(char const *key) const = 0;
        };

        class IDictionary : public IConstDictionary {
        public:
            virtual ~IDictionary() {}

            virtual void setValue(char const *key, char const *value) = 0;
            virtual void clear() = 0;
        };

        class ISequence {
        public:
            virtual ~ISequence() {}

            virtual IConstDictionary const &getMetadata() = 0;
            virtual uint64_t read(char *buffer,
                                  uint64_t buffer_length) = 0;
        };

        class ISequenceIterator {
        public:
            virtual ~ISequenceIterator() {}

            virtual ISequence *nextSequence() = 0;
        };

        class IWriter {
        public:
            virtual ~IWriter() {}

            virtual void createSequence(IConstDictionary const *metadata) = 0;
            virtual void write(char const *buffer,
                               uint64_t length) = 0;
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

#define raise_state(MSG...) {                    \
        raise(INVALID_STATE, MSG);               \
    }

#define implement() {                                                   \
        fprintf(stderr, "implement %s @ %s:%d\n", __FUNCTION__, __FILE__, __LINE__); \
        abort();                                                        \
    }
