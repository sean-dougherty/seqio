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
        seqio_writer_options opts = {SEQIO_FILE_FORMAT_FASTA_GZIP};
        seqio_writer writer;
        seqio_create_writer("/tmp/seqio.fa",
                            opts,
                            &writer);

        seqio_create_sequence(writer);
        seqio_add_metadata(writer, SEQIO_KEY_NAME, "seq1");
        seqio_add_metadata(writer, SEQIO_KEY_COMMENT, "comment1.0 comment1.1");
        seqio_write(writer, seq, seqlen);

        seqio_dispose_writer(&writer);
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


int main(int argc, const char **argv) {
    test_fasta_write();

    test_fasta_plain__sequential();
    test_fasta_gzip__sequential();
    test_fasta_plain__out_of_order();
    test_fasta_gzip__out_of_order();

    cout << "Test successful." << endl;

    return 0;
}
