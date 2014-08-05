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


int main(int argc, const char **argv) {
    test_fasta_plain__sequential();
    test_fasta_gzip__sequential();
    test_fasta_plain__out_of_order();
    test_fasta_gzip__out_of_order();

    cout << "Test successful." << endl;

    return 0;
}
