#include "test_util.hpp"

#include "seqio.h"

using namespace std;


void test_fasta_plain__sequential() {
    seqio_set_err_handler(SEQIO_ERR_HANDLER_ABORT);
    verify_a__sequential("input/a.fa");
}

void test_fasta_gzip__sequential() {
    seqio_set_err_handler(SEQIO_ERR_HANDLER_ABORT);
    verify_a__sequential("input/a.fa.gz");
}

void test_fasta_plain__out_of_order() {
    seqio_set_err_handler(SEQIO_ERR_HANDLER_ABORT);
    verify_a__out_of_order("input/a.fa");
}

void test_fasta_gzip__out_of_order() {
    seqio_set_err_handler(SEQIO_ERR_HANDLER_ABORT);
    verify_a__out_of_order("input/a.fa.gz");
}

void test_fasta_write() {
    seqio_set_err_handler(SEQIO_ERR_HANDLER_ABORT);
    uint32_t const seqlen = 1024 * 1024 * 16;
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
        verify_sequence(sequence, "seq1", "comment1.0 comment1.1", seq);
    }
    
    free(seq);
}

void test_pna_write() {
    seqio_set_err_handler(SEQIO_ERR_HANDLER_ABORT);

    verify_write(SEQIO_FILE_FORMAT_PNA, "/tmp/seqio.pna");
}

void test_read_all__small() {
    seqio_set_err_handler(SEQIO_ERR_HANDLER_ABORT);
    verify_read_all(16);
}

void test_read_all__large() {
    seqio_set_err_handler(SEQIO_ERR_HANDLER_ABORT);
    verify_read_all(16 * 1024 * 1024);
}


int main(int argc, const char **argv) {
    test_pna_write();
    test_fasta_write();

    test_read_all__small();
    test_read_all__large();

    test_fasta_plain__sequential();
    test_fasta_gzip__sequential();

    test_fasta_plain__out_of_order();
    test_fasta_gzip__out_of_order();

    cout << "Test successful." << endl;

    return 0;
}
