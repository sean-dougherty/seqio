#include "fasta.h"
#include "util.h"
#include "kseq.h"

#include <string.h>

#include <iostream>
#include <string>
#include <vector>

// Initialize the kseq library, but disable a warning from it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
KSEQ_INIT(gzFile, gzread)
#pragma GCC diagnostic pop

using namespace std;
using namespace seqio;

void translate(kseq_t *kseq, FastaReader::Translate translate);

void usage(string msg = "") {
    epf("usage: fasta [--translate (none|caps_gatcn)] v[alidate] <fasta...>");
    epf("       fasta [--translate (none|caps_gatcn)] cat <fasta...>");

    if(msg.length() > 0) {
	ep(msg.c_str());
    }

    exit(1);
}

int main(int argc, const char **argv) {

    FastaReader::Translate translate = FastaReader::Translate_None;
    
    int argi = 1;
    for(; argi < argc; argi++) {
	string flag = argv[argi];
	if(flag[0] != '-') break;

	if(flag == "--translate") {
	    argi++;
	    if(argi == argc) usage("Missing translate mode");
	    string transtr = argv[argi];
	    if(transtr == "none") {
		translate = FastaReader::Translate_None;
	    } else if(transtr == "caps_gatcn") {
		translate = FastaReader::Translate_Caps_GATCN;
	    } else {
		usage("Invalid translate mode: "+transtr);
	    }
	} else {
	    usage("Invalid flag: "+flag);
	}
    }

    if((argc - argi) < 1) {
	usage();
    }

    string mode = argv[argi++];
    if((mode == "validate") || (mode == "v")) {
	for(; argi < argc; argi++) {
	    const char *path = argv[argi];

	    cout << "Validating " << path << endl;

	    gzFile fp = gzopen(path, "r");
	    kseq_t *kseq = kseq_init(fp);

	    FastaReader reader(path, translate);
	    FastaSequenceDesc seq;

	    while(true) {
		bool kseqNext = (kseq_read(kseq) >= 0);
		if(kseqNext) {
		    ::translate(kseq, translate);
		}

		bool fastaNext = reader.nextSequence(seq);
		errif(kseqNext != fastaNext, "next mismatch");
		if(!kseqNext) break;

		errif(seq.name != kseq->name.s,
		      "Name mismatch; kseq=%s, fasta=%s.",
		      kseq->name.s, seq.name.c_str());
		errif(seq.comment != kseq->comment.s,
		      "Comment mismatch; kseq=%s, fasta=%s.",
		      kseq->comment.s, seq.comment.c_str());

		cout << "  " << seq.name << " " << seq.comment << endl;

		uint64_t fastaLen = 0;
		uint64_t fastaRc = 0;
		char buf[4 * 1024];
		while( 0 != (fastaRc = reader.read(buf, sizeof(buf))) ) {
		    for(uint64_t i = 0; i < fastaRc; i++) {
			errif(buf[i] != kseq->seq.s[fastaLen + i],
			      "Base mismatch at %ld;"
			      " kseq='%c', fasta='%c'",
			      fastaLen + i,
			      kseq->seq.s[fastaLen + i], buf[i]);
		    }
		    fastaLen += fastaRc;
		}

		errif(fastaLen != kseq->seq.l,
		      "Length mismatch: kseq=%ld fasta=%ld",
		      kseq->seq.l, fastaLen);
	    }

	    kseq_destroy(kseq);
	    gzclose(fp);
	}
	cout << "SUCCESSFULL FASTA VALIDATION." << endl; 
    } else if(mode == "cat") {
	for(; argi < argc; argi++) {
	    const char *path = argv[argi];

	    FastaReader reader(path, translate);
	    FastaSequenceDesc seq;
	    while(reader.nextSequence(seq)) {
		char buf[16 * 1024];
		
		uint64_t rc;
		while( 0 != (rc = reader.read(buf, sizeof(buf))) ) {
		    errif(1 != fwrite(buf, rc, 1, stdout),
			  "Failed writing to stdout");
		}
	    }
	}
    } else {
	usage("Invalid mode: "+mode);
    }
    
    return 0;
}

void translate(kseq_t *kseq, FastaReader::Translate translate) {
    char base_map[256];

    switch(translate) {
    case FastaReader::Translate_None:
	// no-op
	return;
    case FastaReader::Translate_Caps_GATCN:
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
