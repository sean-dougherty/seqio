#include "seqio.h"
#include "util.h"

#include <zlib.h>
#include "kseq.h"

#include <string.h>

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Initialize the kseq library, but disable a warning from it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
KSEQ_INIT(gzFile, gzread)
#pragma GCC diagnostic pop

using namespace std;

void transform(kseq_t *kseq, seqio_base_transform transform);

void usage(string msg = "") {
    epf("usage: fasta [--transform (none|caps_gatcn)] v[alidate] <fasta...>");
    epf("       fasta [--transform (none|caps_gatcn)] cat <fasta...>");
    epf("       fasta [--transform (none|caps_gatcn)] read <fasta...>");
    epf("       fasta [--transform (none|caps_gatcn)] split fasta outdir");

    if(msg.length() > 0) {
        ep(msg.c_str());
    }

    exit(1);
}

int main(int argc, const char **argv) {
    seqio_sequence_options opts = SEQIO_DEFAULT_SEQUENCE_OPTIONS;
    
    int argi = 1;
    for(; argi < argc; argi++) {
        string flag = argv[argi];
        if(flag[0] != '-') break;

        if(flag == "--transform") {
            argi++;
            if(argi == argc) usage("Missing transform mode");
            string transtr = argv[argi];
            if(transtr == "none") {
                opts.base_transform = SEQIO_BASE_TRANSFORM_NONE;
            } else if(transtr == "caps_gatcn") {
                opts.base_transform = SEQIO_BASE_TRANSFORM_CAPS_GATCN;
            } else {
                usage("Invalid transform mode: "+transtr);
            }
        } else {
            usage("Invalid flag: "+flag);
        }
    }

    if((argc - argi) < 1) {
        usage();
    }

    char *buf = nullptr;
    uint64_t buflen;

    string mode = argv[argi++];
    if((mode == "validate") || (mode == "v")) {


        for(; argi < argc; argi++) {
            const char *path = argv[argi];

            cout << "Validating " << path << endl;

            gzFile fp = gzopen(path, "r");
            kseq_t *kseq = kseq_init(fp);

            seqio_sequence_iterator iterator;
            seqio_sequence sequence;
            seqio_create_sequence_iterator(path, opts, &iterator);

            while(true) {
                bool kseqNext = (kseq_read(kseq) >= 0);
                if(kseqNext) {
                    ::transform(kseq, opts.base_transform);
                }

                seqio_next_sequence(iterator, &sequence);
                bool fastaNext = sequence != nullptr;
                errif(kseqNext != fastaNext, "next mismatch");
                if(!kseqNext) break;

                char const *name, *comment;
                seqio_const_dictionary metadata;
                seqio_get_metadata(sequence, &metadata);
                seqio_get_value(metadata, SEQIO_KEY_NAME, &name);
                seqio_get_value(metadata, SEQIO_KEY_COMMENT, &comment);
                errif(string(name) != kseq->name.s,
                      "Name mismatch; kseq=%s, fasta=%s.",
                      kseq->name.s, name);
                errif(string(comment) != kseq->comment.s,
                      "Comment mismatch; kseq=%s, fasta=%s.",
                      kseq->comment.s, comment);

                cout << "  " << name << " " << comment << endl;

                uint64_t seqlen;
                seqio_read_all(sequence, &buf, &buflen, &seqlen);

                errif(seqlen != kseq->seq.l,
                      "Length mismatch: kseq=%zu fasta=%zu",
                      size_t(kseq->seq.l), size_t(seqlen));

                for(uint32_t i = 0; i < seqlen; i++) {
                    errif(buf[i] != kseq->seq.s[i],
                          "Base mismatch at %zu;"
                          " kseq='%c', fasta='%c'",
                          size_t(i),
                          kseq->seq.s[i], buf[i]);
                }

                seqio_dispose_sequence(&sequence);
            }
            kseq_destroy(kseq);
            gzclose(fp);
        }

        cout << "SUCCESSFULL FASTA VALIDATION." << endl; 
    } else if(mode == "cat" || mode == "read") {
        char *buf = nullptr;
        uint64_t buflen;

        for(; argi < argc; argi++) {
            const char *path = argv[argi];

            seqio_sequence_iterator iterator;
            seqio_sequence sequence;
            seqio_create_sequence_iterator(path, opts, &iterator);

            while( (0 == seqio_next_sequence(iterator, &sequence)) && sequence) {
                uint64_t seqlen;
                seqio_read_all(sequence, &buf, &buflen, &seqlen);
                
                if(mode == "cat") {
                    errif(1 != fwrite(buf, seqlen, 1, stdout),
                          "Failed writing to stdout");
                }

                seqio_dispose_sequence(&sequence);
            }

            seqio_dispose_sequence_iterator(&iterator);
        }

    } else if(mode == "split") {
        if( (argc - argi) != 2) {
            usage();
        }
        string path_in = argv[argi++];
        boost::filesystem::path path_outdir = argv[argi++];

        if(!boost::filesystem::exists(path_outdir)) {
            errif(!boost::filesystem::create_directories(path_outdir), "Failed creating output directory.");
        }
        errif(!boost::filesystem::is_directory(path_outdir), "Output not a directory");

        seqio_sequence_iterator iterator;
        seqio_sequence sequence;
        seqio_create_sequence_iterator(path_in.c_str(), opts, &iterator);

        while( (0 == seqio_next_sequence(iterator, &sequence)) && sequence) {
            char const *name, *comment;
            seqio_const_dictionary metadata;
            seqio_get_metadata(sequence, &metadata);
            seqio_get_value(metadata, SEQIO_KEY_NAME, &name);
            seqio_get_value(metadata, SEQIO_KEY_COMMENT, &comment);

            uint64_t seqlen;
            seqio_read_all(sequence, &buf, &buflen, &seqlen);

            boost::filesystem::path path_out = path_outdir / (string(name) + ".fa");

            seqio_writer writer;
            seqio_create_writer(path_out.c_str(),
                                SEQIO_DEFAULT_WRITER_OPTIONS,
                                &writer);

            seqio_create_sequence(writer, metadata);
            seqio_write(writer, buf, seqlen);

            seqio_dispose_writer(&writer);
        }
    } else {
        usage("Invalid mode: "+mode);
    }

    seqio_dispose_buffer(&buf);
    
    return 0;
}

void transform(kseq_t *kseq, seqio_base_transform base_transform) {
    char base_map[256];

    switch(base_transform) {
    case SEQIO_BASE_TRANSFORM_NONE:
        // no-op
        return;
    case SEQIO_BASE_TRANSFORM_CAPS_GATCN:
        for(int i = 0; i < 256; i++) {
            base_map[i] = 'N';
        }
        base_map['A'] = 'A';
        base_map['a'] = 'A';
        base_map['C'] = 'C';
        base_map['c'] = 'C';
        base_map['G'] = 'G';
        base_map['g'] = 'G';
        base_map['T'] = 'T';
        base_map['t'] = 'T';
        break;
    default:
        err("Tranlate mode not implemented for kseq");
        break;
    }

    for(size_t i = 0; i < kseq->seq.l; i++) {
        kseq->seq.s[i] = base_map[(uint8_t)kseq->seq.s[i]];
    }
}
