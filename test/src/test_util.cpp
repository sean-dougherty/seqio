#include "test_util.hpp"

#include <assert.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/resource.h>

#include <cstring>
#include <random>
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

char *create_random_bases(uint64_t len, int seed) {
    char *seq = (char *)malloc(len + 1);
    char bases[] = {'A','T','G','C'};
    
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> dist(0, 3);

    for(uint64_t i = 0; i < len; i++) {
        seq[i] = bases[dist(generator)];
    }
    seq[len] = 0;

    return seq;
}

void write_file(char const *path, vector<seqspec_t> seqs) {
    seqio_writer writer;
    seqio_create_writer(path,
                        SEQIO_DEFAULT_WRITER_OPTIONS,
                        &writer);

    seqio_dictionary metadata;
    seqio_create_dictionary(&metadata);

    for(seqspec_t spec: seqs) {
        seqio_set_value(metadata, SEQIO_KEY_NAME, spec.name);
        seqio_set_value(metadata, SEQIO_KEY_COMMENT, spec.comment);

        seqio_create_sequence(writer, metadata);
        seqio_write(writer, spec.bases, strlen(spec.bases));
    }

    seqio_dispose_writer(&writer);
}

void verify_basic_metadata(seqio_sequence sequence, char const *name, char const *comment) {
    seqio_const_dictionary dict;
    seqio_get_metadata(sequence, &dict);

    {
        const char *name_;
        seqio_get_value(dict, SEQIO_KEY_NAME, &name_);
        assert(0 == strcmp(name, name_));
    }

    {    
        const char *comment_;
        seqio_get_value(dict, SEQIO_KEY_COMMENT, &comment_);
        assert(0 == strcmp(comment_, comment));
    }
}

void verify_bases(seqio_sequence sequence, char const *expected, uint64_t buflen) {
    uint64_t seqlen = strlen(expected);
    char *buf = (char *)malloc(seqlen+1);
    memset(buf, 0, seqlen+1);

    uint64_t total_read_length = 0;
    uint64_t read_length;
    while( (seqio_read(sequence, buf + total_read_length, buflen, &read_length) == SEQIO_SUCCESS)
           && (read_length > 0) ) {

        for(uint64_t i = total_read_length; i < total_read_length + read_length; i++) {
            assert(buf[i] == expected[i]);
        }

        total_read_length += read_length;
        assert(total_read_length <= seqlen);
    }

    assert(total_read_length == seqlen);
    assert( 0 == strcmp(expected, buf) );

    seqio_read(sequence, buf, buflen, &read_length);
    assert(read_length == 0);

    free(buf);
}

void verify_sequence(char const *path,
                     seqio_base_transform base_transform,
                     char const *name,
                     char const *comment,
                     char const *bases) {
    seqio_sequence_options options = SEQIO_DEFAULT_SEQUENCE_OPTIONS;
    options.base_transform = base_transform;

    seqio_sequence_iterator iterator;
    seqio_create_sequence_iterator(path,
                                   options,
                                   &iterator);

    {
        seqio_sequence sequence;
        seqio_next_sequence(iterator, &sequence);
        verify_sequence(sequence, name, comment, bases);
    }
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

void verify_read_all(uint64_t seqlen) {
    char *seq = create_random_bases(seqlen);

    {
        seqio_writer_options opts = {SEQIO_FILE_FORMAT_FASTA};
        seqio_writer writer;
        seqio_create_writer("/tmp/seqio.fa",
                            opts,
                            &writer);

        seqio_dictionary metadata;
        seqio_create_dictionary(&metadata);

        seqio_set_value(metadata, SEQIO_KEY_NAME, "seq1");
        seqio_set_value(metadata, SEQIO_KEY_COMMENT, "comment1.0 comment1.1");
        seqio_create_sequence(writer, metadata);
        seqio_write(writer, seq, seqlen);

        seqio_dispose_writer(&writer);
        seqio_dispose_dictionary(&metadata);
    }

    {
        seqio_sequence_iterator iterator;
        seqio_create_sequence_iterator("/tmp/seqio.fa",
                                       SEQIO_DEFAULT_SEQUENCE_OPTIONS,
                                       &iterator);

        seqio_sequence sequence;
        seqio_next_sequence(iterator, &sequence);
        
        char *buffer = nullptr;
        uint64_t buffer_length;
        uint64_t read_length;
        
        seqio_read_all(sequence, &buffer, &buffer_length, &read_length);

        assert(read_length == seqlen);
        assert(buffer_length >= seqlen + 1);
        assert(strlen(buffer) == seqlen);
        assert(0 == strcmp(buffer, seq));

        seqio_dispose_buffer(&buffer);
    }
    
    free(seq);
}

void verify_write(seqio_file_format file_format, char const *path) {
    uint64_t const seqlen = 1024 * 1024 * 16;
    char *seq = create_random_bases(seqlen);

    {
        seqio_writer_options opts = {file_format};
        seqio_writer writer;
        seqio_create_writer(path,
                            opts,
                            &writer);

        seqio_dictionary metadata;
        seqio_create_dictionary(&metadata);

        seqio_set_value(metadata, SEQIO_KEY_NAME, "seq1");
        seqio_set_value(metadata, SEQIO_KEY_COMMENT, "comment1.0 comment1.1");
        seqio_create_sequence(writer, metadata);
        seqio_write(writer, seq, seqlen);

        seqio_dispose_writer(&writer);
    }

    {
        seqio_sequence_iterator iterator;
        seqio_create_sequence_iterator(path,
                                       SEQIO_DEFAULT_SEQUENCE_OPTIONS,
                                       &iterator);

        seqio_sequence sequence;
        seqio_next_sequence(iterator, &sequence);
        verify_sequence(sequence, "seq1", "comment1.0 comment1.1", seq);
    }

    free(seq);    
}
