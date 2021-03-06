#include "seqio.h"

#include "fasta.hpp"
#include "pna_impl.hpp"

#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace seqio::impl;
using std::map;
using std::string;

seqio_sequence_options const SEQIO_DEFAULT_SEQUENCE_OPTIONS = {
    SEQIO_FILE_FORMAT_DEDUCE,
    SEQIO_BASE_TRANSFORM_NONE
};

seqio_writer_options const SEQIO_DEFAULT_WRITER_OPTIONS = {
    SEQIO_FILE_FORMAT_DEDUCE
};

static seqio_status  __err_abort(seqio_err_info err_info);
static seqio_status  __err_exit(seqio_err_info err_info);
static seqio_status __err_return(seqio_err_info err_info);
seqio_err_handler const SEQIO_ERR_HANDLER_ABORT = &__err_abort;
seqio_err_handler const SEQIO_ERR_HANDLER_EXIT = &__err_exit;
seqio_err_handler const SEQIO_ERR_HANDLER_RETURN = &__err_return;
static seqio_err_handler err_handler = SEQIO_ERR_HANDLER_EXIT;

#define INITIAL_READ_ALL_BUFFER_LENGTH (16 * 1024)

#define err(STATUS, MSG...) {                                       \
        char message[4096];                                         \
        sprintf(message, MSG);                                      \
        seqio_err_info err_info = {SEQIO_ERR_##STATUS, message};    \
        return err_handler(err_info);                               \
    }

#define err_parm(MSG) {                                     \
        err(INVALID_PARAMETER, MSG);                        \
    }

#define err_fnf(MSG...) {                                   \
        err(FILE_NOT_FOUND, MSG);                           \
    }

#define check_null(PARAMETER) {                             \
        if((PARAMETER) == nullptr) {                        \
            err_parm(#PARAMETER " is null");                \
        }                                                   \
    }

#define TRY(STMT,CATCH_STMT)                    \
    try {                                       \
    STMT;                                       \
    } catch(seqio::impl::Exception x) {         \
        CATCH_STMT;                             \
    return err_handler(x.err_info);             \
    }

seqio_status seqio_create_sequence_iterator(char const *path,
                                            seqio_sequence_options options,
                                            seqio_sequence_iterator *iterator) {
    check_null(path);
    check_null(iterator);

    {
        struct stat buf;
        if(0 != stat(path, &buf)) {
            err_fnf("No such file: %s", path);
        }
    }

    ISequenceIterator *impl;

    if(options.file_format == SEQIO_FILE_FORMAT_DEDUCE) {
        if(is_pna_file_content(path)) {
            options.file_format = SEQIO_FILE_FORMAT_PNA;
        } else {
            options.file_format = SEQIO_FILE_FORMAT_FASTA;
        }
    }

    try {
        switch(options.file_format) {
        case SEQIO_FILE_FORMAT_FASTA:
        case SEQIO_FILE_FORMAT_FASTA_GZIP:
            impl = new FastaSequenceIterator(path, options.base_transform);
            break;
        case SEQIO_FILE_FORMAT_PNA:
            impl = new PnaSequenceIterator(path, options.base_transform);
            break;
        default:
            raise_parm("Invalid file format.");
        }
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    *iterator = (seqio_sequence_iterator)impl;

    return SEQIO_SUCCESS;
}

seqio_status seqio_dispose_sequence_iterator(seqio_sequence_iterator *iterator) {
    if(iterator && *iterator) {
        try {
            delete (ISequenceIterator *)*iterator;
        } catch(Exception x) {
            return err_handler(x.err_info);
        }
        *iterator = nullptr;
    }

    return SEQIO_SUCCESS;
}

seqio_status seqio_next_sequence(seqio_sequence_iterator iterator,
                                 seqio_sequence *sequence) {
    check_null(iterator);
    check_null(sequence);

    try {
        *sequence = (seqio_sequence)((ISequenceIterator *)iterator)->nextSequence();
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    return SEQIO_SUCCESS;    
}

seqio_status seqio_dispose_sequence(seqio_sequence *sequence) {
    if(sequence && *sequence) {
        try {
            delete (ISequence *)*sequence;
        } catch(Exception x) {
            return err_handler(x.err_info);
        }
        *sequence = nullptr;
    }

    return SEQIO_SUCCESS;
}

seqio_status seqio_get_metadata(seqio_sequence sequence,
                                seqio_const_dictionary *dict) {
    check_null(sequence);
    check_null(dict);

    try {
        ISequence *iseq = (ISequence *)sequence;
        *dict = (seqio_const_dictionary)&iseq->getMetadata();
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    return SEQIO_SUCCESS;
}

seqio_status seqio_read(seqio_sequence sequence,
                        char *buffer,
                        uint64_t buffer_length,
                        uint64_t *read_length) {
    check_null(sequence);
    check_null(buffer);
    check_null(read_length);

    try {
        uint64_t n = ((ISequence *)sequence)->read(buffer, buffer_length);
        *read_length = n;
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    return SEQIO_SUCCESS;    
}

seqio_status seqio_read_all(seqio_sequence sequence,
                            char **buffer,
                            uint64_t *buffer_length,
                            uint64_t *read_length) {
    check_null(sequence);
    check_null(buffer);
    check_null(buffer_length);
    check_null(read_length);

    ISequence *iseq = (ISequence *)sequence;

    try {
        if(*buffer == nullptr) {
            *buffer_length = INITIAL_READ_ALL_BUFFER_LENGTH;
            *buffer = (char *)malloc(*buffer_length);
            if(*buffer == nullptr)
                raise_oom("Cannot allocate %zu bytes", size_t(*buffer_length));
        }

        auto grow_to = [=] (uint64_t new_length) {
            char *buffer_ = (char *)realloc(*buffer, new_length);
            if(buffer_ == nullptr) {
                raise_oom("Cannot allocate %zu bytes", size_t(new_length));
            }
            *buffer = buffer_;
            *buffer_length = new_length;
        };

        uint64_t nread;
        uint64_t nread_total = 0;

        while( (nread = iseq->read(*buffer + nread_total, *buffer_length - nread_total)) != 0) {
            nread_total += nread;
            if(nread_total == *buffer_length) {
                grow_to(*buffer_length * 2);
            }
        }

        if(nread_total == *buffer_length) {
            grow_to(*buffer_length + 1);
        }
        (*buffer)[nread_total] = '\0';

        *read_length = nread_total;
    } catch(Exception x) {
        *read_length = 0;
        return err_handler(x.err_info);
    }
    
    return SEQIO_SUCCESS;    
}

seqio_status seqio_create_writer(char const *path,
                                 seqio_writer_options options,
                                 seqio_writer *writer) {
    check_null(path);
    check_null(writer);

    IWriter *iwriter;

    try {
        if(options.file_format == SEQIO_FILE_FORMAT_DEDUCE) {
            if(is_fasta_file_name(path))
                options.file_format = SEQIO_FILE_FORMAT_FASTA;
            else if(is_fasta_gzip_file_name(path))
                options.file_format = SEQIO_FILE_FORMAT_FASTA_GZIP;
            else if(is_pna_file_name(path))
                options.file_format = SEQIO_FILE_FORMAT_PNA;
            else
                raise_parm("Cannot deduce file format from file extension.");
        }

        switch(options.file_format) {
        case SEQIO_FILE_FORMAT_FASTA:
        case SEQIO_FILE_FORMAT_FASTA_GZIP:
            iwriter = new FastaWriter(path, options.file_format);
            break;
        case SEQIO_FILE_FORMAT_PNA:
            iwriter = new PnaWriter(path, options.file_format);
            break;
        default:
            raise_parm("Invalid file_format specified.");
        }

    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    *writer = (seqio_writer)iwriter;

    return SEQIO_SUCCESS;
}

seqio_status seqio_dispose_writer(seqio_writer *writer) {
    if(writer && *writer) {
        try {
            delete (IWriter *)*writer;
            *writer = nullptr;
        } catch(Exception x) {
            return err_handler(x.err_info);
        }
    }

    return SEQIO_SUCCESS;
}

seqio_status seqio_create_sequence(seqio_writer writer,
                                   seqio_const_dictionary dict) {
    check_null(writer);
    check_null(dict);

    try {
        ((IWriter *)writer)->createSequence((IConstDictionary const *)dict);
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    return SEQIO_SUCCESS;
}

seqio_status seqio_write(seqio_writer writer,
                         char const *buffer,
                         uint64_t length) {
   check_null(writer);
   check_null(buffer);

    try {
        ((IWriter *)writer)->write(buffer, length);
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    return SEQIO_SUCCESS;
}

class Dictionary : public IDictionary {
public:
    virtual ~Dictionary() {}

    virtual bool hasKey(char const *key) const override {
        return dict.find(key) != dict.end();
    }

    virtual uint32_t getKeyCount() const override {
        return dict.size();
    }

    virtual char const *getKey(uint32_t key_index) const override {
        if(key_index >= dict.size())
            raise_parm("Index exceeds dictionary size.");

        auto it = dict.begin();
        for(uint32_t i = 0; i < key_index; i++)
            ++it;

        return it->first.c_str();
    }

    virtual char const *getValue(char const *key) const override {
        auto it = dict.find(key);
        if(it == dict.end())
            raise_parm("Unknown key: %s", key);

        return it->second.c_str();
    }

    virtual void setValue(char const *key, char const *value) override {
        dict[key] = value;
    }

    virtual void clear() override {
        dict.clear();
    }

private:
    map<string, string> dict;
};

seqio_status seqio_create_dictionary(seqio_dictionary *dict) {
   check_null(dict);

   *dict = (seqio_dictionary)new Dictionary();

    return SEQIO_SUCCESS;
}

seqio_status seqio_dispose_dictionary(seqio_dictionary *dict) {
    if(dict && *dict) {
        delete (Dictionary *)*dict;
        *dict = nullptr;
    }

    return SEQIO_SUCCESS;
}

seqio_status seqio_set_value(seqio_dictionary dict,
                             char const *key,
                             char const *value) {
    check_null(dict);
    check_null(key);
    check_null(value);
    
    try {
        ((IDictionary *)dict)->setValue(key, value);
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    return SEQIO_SUCCESS;        
}

seqio_status seqio_clear(seqio_dictionary dict) {
    check_null(dict);

    try {
        ((IDictionary *)dict)->clear();
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    return SEQIO_SUCCESS;        
}

seqio_status seqio_get_key_count(seqio_const_dictionary dict,
                                 uint32_t *count) {
    check_null(dict);
    check_null(count);
    
    try {
        *count = ((IConstDictionary const *)dict)->getKeyCount();
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    return SEQIO_SUCCESS;        
}

seqio_status seqio_get_key(seqio_const_dictionary dict,
                           uint32_t key_index,
                           char const **key) {
    check_null(dict);
    check_null(key);
    
    try {
        *key = ((IConstDictionary const *)dict)->getKey(key_index);
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    return SEQIO_SUCCESS;        
}


seqio_status seqio_get_value(seqio_const_dictionary dict,
                             char const *key,
                             char const **value) {
    check_null(dict);
    check_null(key);
    check_null(value);
    
    try {
        *value = ((IConstDictionary const *)dict)->getValue(key);
    } catch(Exception x) {
        return err_handler(x.err_info);
    }

    return SEQIO_SUCCESS;        
}

seqio_status seqio_dispose_buffer(char **buffer) {
    if(buffer && *buffer) {
        free(*buffer);
        *buffer = NULL;
    }

    return SEQIO_SUCCESS;
}

seqio_status seqio_set_err_handler(seqio_err_handler err_handler_) {
    check_null(err_handler_);

    err_handler = err_handler_;

    return SEQIO_SUCCESS;
}

static seqio_status __err_abort(seqio_err_info err_info) {
    fprintf(stderr, "seqio error: status=%d, message=%s\n",
            (int)err_info.status,
            err_info.message);
    abort();
}

static seqio_status __err_exit(seqio_err_info err_info) {
    fprintf(stderr, "seqio error: status=%d, message=%s\n",
            (int)err_info.status,
            err_info.message);
    exit(err_info.status);
}

static seqio_status __err_return(seqio_err_info err_info) {
    return err_info.status;
}
