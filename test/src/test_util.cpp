#include "test_util.hpp"

#include <assert.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/resource.h>

#include <cstring>
#include <vector>

using namespace std;


int open_count(char const *path) {
    int i = 0;
    int fd_counter = 0;
    int max_fd_number = 0;
    struct stat   stat_path, stat_fd;
    struct rlimit rlimits;

    assert( 0 == stat(path, &stat_path) );
 
    max_fd_number = getdtablesize();
    getrlimit(RLIMIT_NOFILE, &rlimits);
 
    for ( i = 0; i <= max_fd_number; i++ ) {
        if( (0 == fstat(i, &stat_fd)) && (stat_fd.st_ino == stat_path.st_ino) ) {
            fd_counter++;
        }
    }

    return fd_counter;
}

void verify_basic_metadata(seqio_sequence sequence, char const *name, char const *comment) {
    {    
        const char *name_;
        seqio_get_value(sequence, SEQIO_KEY_NAME, &name_);
        assert(0 == strcmp(name, name_));
    }

    {    
        const char *comment_;
        seqio_get_value(sequence, SEQIO_KEY_COMMENT, &comment_);
        assert(0 == strcmp(comment_, comment));
    }
}

void verify_bases(seqio_sequence sequence, char const *expected, uint32_t buflen) {
    uint32_t seqlen = strlen(expected);
    char buf[seqlen+1];
    memset(buf, 0, seqlen+1);

    uint32_t total_read_length = 0;
    uint32_t read_length;
    while( (seqio_read(sequence, buf + total_read_length, buflen, &read_length) == SEQIO_SUCCESS)
           && (read_length > 0) ) {

        total_read_length += read_length;
        assert(total_read_length <= seqlen);
    }

    assert(total_read_length == seqlen);
    assert( 0 == strcmp(expected, buf) );

    seqio_read(sequence, buf, buflen, &read_length);
    assert(read_length == 0);
}

void verify_sequence(seqio_sequence sequence,
                     char const *name,
                     char const *comment,
                     char const *bases) {
    verify_basic_metadata(sequence, name, comment);
    verify_bases(sequence, bases, 2);

    seqio_dispose_sequence(&sequence);
    assert(!sequence);
}

void verify_a__sequential(char const *path) {
    seqio_sequence_iterator iterator;
    seqio_create_sequence_iterator(path,
                                   SEQIO_DEFAULT_SEQUENCE_OPTIONS,
                                   &iterator);

    {
        seqio_sequence sequence;
        seqio_next_sequence(iterator, &sequence);
        verify_sequence(sequence, "seq1", "comment1.0 comment1.1", "aAgGcCtT");
    }

    {
        seqio_sequence sequence;
        seqio_next_sequence(iterator, &sequence);
        verify_sequence(sequence, "seq2", "", "acgtACGT");
    }
 
    seqio_dispose_sequence_iterator(&iterator);
    assert(!iterator);

    assert( 0 == open_count(path) );
}

void verify_a__out_of_order(char const *path) {
    assert( 0 == open_count(path) );

    seqio_sequence_iterator iterator;
    seqio_create_sequence_iterator(path,
                                   SEQIO_DEFAULT_SEQUENCE_OPTIONS,
                                   &iterator);

    vector<seqio_sequence> sequences;
    {
        seqio_sequence sequence;
        seqio_next_sequence(iterator, &sequence);
        sequences.push_back(sequence);

        seqio_next_sequence(iterator, &sequence);
        sequences.push_back(sequence);
    }
    seqio_dispose_sequence_iterator(&iterator);

    verify_sequence(sequences[1], "seq2", "", "acgtACGT");
    verify_sequence(sequences[0], "seq1", "comment1.0 comment1.1", "aAgGcCtT");
}
