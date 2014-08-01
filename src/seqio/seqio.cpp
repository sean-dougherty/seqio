#include "seqio.h"

#include <stdio.h>
#include <stdlib.h>

static seqio_status  __err_abort(seqio_err_info err_info);
static seqio_status __err_return(seqio_err_info err_info);
seqio_err_handler const SEQIO_ERR_HANDLER_ABORT = &__err_abort;
seqio_err_handler const SEQIO_ERR_HANDLER_RETURN = &__err_return;
static seqio_err_handler err_handler = SEQIO_ERR_HANDLER_ABORT;

enum class Format { PNA, FASTA };

struct __sequence_iterator {
    Format format;
    void *impl;

    __sequence_iterator(Format format_, void *impl_)
        : format(format_), impl(impl_) {
    }
};

#define err(STATUS, MSG...) {                                       \
        char message[4096];                                         \
        sprintf(message, MSG);                                      \
        seqio_err_info err_info = {SEQIO_ERR_##STATUS, message};    \
        return err_handler(err_info);                               \
    }

#define bad_parm(MSG) {                                     \
        err(INVALID_PARAMETER, MSG);                        \
    }

#define check_null(PARAMETER) {                             \
        if((PARAMETER) == nullptr) {                        \
            bad_parm(#PARAMETER " is null");                \
        }                                                   \
    }

extern "C" {

    seqio_status seqio_get_sequences(char const *path,
                                     seqio_sequence_options options,
                                     seqio_sequence_iterator *iterator) {
        check_null(path);
        check_null(iterator);

        *iterator = new __sequence_iterator(Format::FASTA, nullptr);

        return SEQIO_SUCCESS;
    }

}

static seqio_status __err_abort(seqio_err_info err_info) {
    fprintf(stderr, "seqio error: status=%d, message=%s\n",
            (int)err_info.status,
            err_info.message);
    exit(1);
}

static seqio_status __err_return(seqio_err_info err_info) {
    return err_info.status;
}
