#include "fasta.h"
#include "pna.h"
#include "util.h"

#include <string.h>

#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace seqio;

shared_ptr<PnaWriter> create_writer(const string &path);
void write_seq(FastaReader &reader,
               FastaSequenceDesc &seq, // todo: keep "current seq" in reader
               const string &path_fasta, // todo: put this and index in the desc
               int index_fasta,
               shared_ptr<PnaWriter> fwriter);
void create_file_metadata(shared_ptr<PnaWriter> fwriter);
void create_seq_metadata(shared_ptr<PnaSequenceWriter> writer,
                         const string &path_fasta,
                         int index_fasta,
                         FastaSequenceDesc &seq);
bool parse_ncbi_name(const string &name, map<string, string> &result);
bool parse_ncbi_comment(const string &comment, map<string, string> &result);

void usage(string msg = "") {
    epf("usage: pna a[ssemble] <output_dir> <input_fasta...>");
    epf("       pna c[reate] [-q] <output_pna> <input_fasta...>");
    epf("       pna t[able] [-s] <pna...>");
    epf("       pna v[alidate] [--seek --buflen len] <pna...>");
    epf("       pna cat <pna...>");
    epf("       pna f[asta] [--noN] <pna...>");
    epf("       pna r[ead] [--packed] <pna...>");

    if(msg.length() > 0) {
        ep(msg.c_str());
    }

    exit(1);
}

int main(int argc, const char **argv) {
    if(argc < 2) {
        usage();
    }

    int argi = 1;
    string mode = argv[argi++];

    if((mode == "a") || (mode == "assemble")) {
        if((argc - argi) < 2) {
            usage("Missing assemble arguments");
        }

        map<string, shared_ptr<PnaWriter> > fwriters;
        
        string output_dir = argv[argi++];
        for(; argi < argc; argi++) {
            string path_fasta = argv[argi];

            cout << "Processing " << path_fasta << "..." << endl;
            
            FastaReader reader(path_fasta.c_str(),
                               FastaReader::Translate_Caps_GATCN);
            FastaSequenceDesc seq;
            for(int iseq = 0; reader.nextSequence(seq); iseq++) {
                map<string,string> attrs;
                errif(!parse_ncbi_comment(seq.comment, attrs),
                      "Cannot determine assembly of %s in %s\n"
                      "comment=%s",
                      seq.name.c_str(), path_fasta.c_str(), seq.comment.c_str());

                string assembly = join(split(attrs["assembly"]), "-");
                string species = join(split(attrs["species"]), "_");
                string path_assembly = pathcat(output_dir, species+"_"+assembly+".pna");
                shared_ptr<PnaWriter> fwriter = fwriters[path_assembly];
                if(!fwriter) {
                    fwriter = create_writer(path_assembly);
                    fwriters[path_assembly] = fwriter;
                }

                write_seq(reader, seq, path_fasta, iseq, fwriter);
            }
        }
    } else if((mode == "c") || (mode == "create")) {
        if((argc - argi) < 2) {
            usage("Missing create arguments");
        }
        bool quiet = false;
        if(0 == strcmp(argv[argi], "-q")) {
            quiet = true;
            argi++;
        }
        string path_pna = argv[argi++];
        shared_ptr<PnaWriter> fwriter = create_writer(path_pna);

        for(; argi < argc; argi++) {
            string path_fasta = argv[argi];
            if(!quiet) cout << "Importing " << path_fasta << endl;

            FastaReader reader(path_fasta.c_str(),
                               FastaReader::Translate_Caps_GATCN);
            FastaSequenceDesc seq;
            for(int iseq = 0; reader.nextSequence(seq); iseq++) {
                if(!quiet) cout << "  " << seq.name << " " << seq.comment << endl;
                write_seq(reader, seq, path_fasta, iseq, fwriter);
            }
        }
    } else if((mode == "t") || (mode == "table")) {
        bool summary = false;
        for(; argi < argc; argi++) {
            string arg = argv[argi];
            if(arg[0] != '-')
                break;

            if(arg == "-s") {
                summary = true;
            } else {
                usage("Invalid table flag: " + arg);
            }
        }

        for(; argi < argc; argi++) {
            const char *path = argv[argi];

            cout << path << endl;

            PnaReader reader(path);

            if(!summary) {
                for(uint64_t i = 0; i < reader.getSequenceCount(); i++) {
                    shared_ptr<PnaSequenceReader> seq =
                        reader.openSequence(i);

                    cout << "  sequence " << (i+1) << endl;
                    const PnaMetadata metadata = seq->getMetadata();
                    const char *key, *value;
                    for(uint32_t i = 0; metadata.pair(i, &key, &value); i++) {
                        cout << "    '" << key << "' --> '" << value << "'" << endl;
                    }
                }
            }

            cout << "  Sequence count: " << reader.getSequenceCount() << endl;
            cout << "  Max fragments: " << reader.getMaxSequenceFragments() << endl;
            cout << "  Max packed bytes: " << reader.getMaxPackedBasesLength() << endl;
        }
    } else if((mode == "validate") || (mode == "v")) {
        bool seek = false;
        uint64_t buflen = 0;
        for(; argi < argc; argi++) {
            string arg = argv[argi];
            if(arg[0] != '-')
                break;

            if(arg == "--seek") {
                cerr << "!!!" << endl;
                cerr << "!!! WARNING! SEEK VALIDATION WILL TAKE A LONG TIME!" << endl;
                cerr << "!!!" << endl;
                seek = true;
            } else if(arg == "--buflen") {
                if(++argi == argc) usage("Missing --buflen arg");
                buflen = uint64_t(atol(argv[argi]));
                cout << "Using buflen = " << buflen << endl;
            } else {
                usage("Invalid validate flag: " + arg);
            }
        }

        for(; argi < argc; argi++) {
            const char *path_pna = argv[argi];
            cout << "Validating " << path_pna << endl;

            PnaReader pnaReader(path_pna);
            shared_ptr<FastaReader> fastaReader;
            for(uint64_t i = 0; i < pnaReader.getSequenceCount(); i++) {
                PnaMetadata metadata = pnaReader.getSequenceMetadata(i);
                const char *format = metadata.value("origin.format");
                errif(0 != strcmp(format, "fasta"),
                      "Unimplement validate format: %s", format);
                const char *path_fasta = metadata.value("origin.path");
                int index = atoi(metadata.value("origin.index"));
                errif(index < 0, "Have more than 2GB sequences in fasta?!");
                const char *name = metadata.value("fasta.name");
                const char *comment = metadata.value("fasta.comment");

                if(index == 0) {
                    cout << "  " << path_fasta << endl;
                    fastaReader.reset(new FastaReader(path_fasta,
                                                      FastaReader::Translate_Caps_GATCN));
                }

                FastaSequenceDesc fastaSeq;
                errif(!fastaReader->nextSequence(fastaSeq),
                      "No such fasta sequence: %s:%d",
                      path_fasta, index);

                cout << "    " << name << " " << comment << endl;

                errif(fastaSeq.name != name,
                      "Name mismatch: fasta=%s, pna=%s",
                      fastaSeq.name.c_str(), name);
                errif(fastaSeq.comment != comment,
                      "Comment mismatch: fasta=%s, pna=%s",
                      fastaSeq.comment.c_str(), comment);

                shared_ptr<PnaSequenceReader> pnaSeq =
                    pnaReader.openSequence(i);
                
                if(!seek) {
                    if(buflen == 0)
                        buflen = 4*1024;

                    uint64_t fastaLen = 0;
                    uint64_t pnaLen = 0;
                    uint64_t fastaRc;
                    uint64_t pnaRc;
                    char fastaBuf[buflen];
                    char pnaBuf[buflen];

                    while(true) {
                        fastaRc = fastaReader->read(fastaBuf, sizeof(fastaBuf));
                        pnaRc = pnaSeq->read(pnaBuf, sizeof(pnaBuf));
                        errif(fastaRc != pnaRc,
                              "Read rc mismatch for %s:%s at offset %ld;"
                              "fasta=%ld, pna=%ld",
                              path_fasta, name, fastaLen,
                              fastaRc, pnaRc);

                        if(fastaRc == 0) break;

                        for(uint64_t j = 0; j < fastaRc; j++) {
                            errif(fastaBuf[j] != pnaBuf[j],
                                  "Base mismatch for %s:%s at offset %ld;"
                                  "fasta='%c', pna='%c'",
                                  path_fasta, name, fastaLen + j,
                                  fastaBuf[j], pnaBuf[j]);
                        }

                        fastaLen += fastaRc;
                        pnaLen += pnaRc;
                    }

                    errif(fastaLen != pnaSeq->size(),
                          "Incorrect sequence length reported by PnaSequenceReader;"
                          "expected=%ld, actual=%ld",
                          fastaLen, pnaSeq->size());
                } else {
                    if(buflen == 0)
                        buflen = 16;
                    uint64_t seqOff = 0;

                    cout << "  sequence length=" << pnaSeq->size() << endl;

                    for(int it = 0; ; it++) {
                        if((it % 1000) == 0) {
                            cout << "  offset=" << seqOff << endl;
                        }

                        char fastaBuf[4*1024];
                        uint64_t fastaRc = fastaReader->read(fastaBuf, sizeof(fastaBuf));
                        if(fastaRc == 0)
                            break;

                        for(uint64_t off = fastaRc - 1; ; off--) {
                            uint64_t nread = min(buflen, fastaRc - off);
                            char pnaBuf[nread];

                            pnaSeq->seek(seqOff + off);
                            uint64_t pnaRc = pnaSeq->read(pnaBuf, nread);
                            errif(pnaRc != nread,
                                  "Incomplete read at %lu:%lu; "
                                  "expecting %ld, found %ld",
                                  seqOff, off, nread, pnaRc);

                            for(uint64_t j = 0; j < nread; j++) {
                                errif(fastaBuf[off + j] != pnaBuf[j],
                                      "Base at offset %lu:%lu+%lu; "
                                      "fasta='%c', pna='%c'\n"
                                      "fasta: %s\n"
                                      "pna:   %s\n",
                                      seqOff, off, j,
                                      fastaBuf[off + j], pnaBuf[j],
                                      strndup(fastaBuf + off, nread),
                                      strndup(pnaBuf, nread));
                            }

                            if(off == 0)
                                break;
                        }
                        seqOff += fastaRc;
                    }
                }
            }

            cout << "SUCCESSFUL VALIDATION OF " << path_pna << endl;
        }
    } else if(mode == "cat") {
        bool packed = false;
        for(; argi < argc; argi++) {
            string arg = argv[argi];
            if(arg[0] != '-')
                break;

            if(arg == "--packed") {
                packed = true;
            } else {
                usage("Invalid validate flag: " + arg);
            }
        }
        for(; argi < argc; argi++) {
            const char *path = argv[argi];

            PnaReader reader(path);
            if(!packed) {
                for(uint64_t i  = 0; i < reader.getSequenceCount(); i++) {
                    shared_ptr<PnaSequenceReader> seq = reader.openSequence(i);
                    char buf[16 * 1024];
                
                    uint64_t rc;
                    while( 0 != (rc = seq->read(buf, sizeof(buf))) ) {
                        errif(1 != fwrite(buf, rc, 1, stdout),
                              "Failed writing to stdout");
                    }
                }
            } else {
                uint64_t buflen = reader.getMaxPackedBasesLength();
                uint8_t *buf = new uint8_t[buflen];
                PnaSequenceReader::packed_read_result_t result;

                for(uint64_t i  = 0; i < reader.getSequenceCount(); i++) {
                    shared_ptr<PnaSequenceReader> seq = reader.openSequence(i);
                    result = seq->packed_read(buf, buflen);
                    errif(1 != fwrite(result.seqfragments.begin,
                                      sizeof(seqfragment_t) * result.seqfragments.count,
                                      1,
                                      stdout),
                          "Failed writing seqfragments");
                    errif(1 != fwrite(result.packed_bases.begin,
                                      result.packed_bases.buflen,
                                      1,
                                      stdout),
                          "Failed writing seqfragments");
                }
            }
        }
    } else if((mode == "fasta") || (mode == "f")) {
        uint32_t flags = PnaSequenceReader::Standard;
        for(; argi < argc; argi++) {
            string arg = argv[argi];
            if(arg[0] != '-')
                break;

            if(arg == "--noN") {
                flags |= PnaSequenceReader::IgnoreN;
            } else {
                usage("Invalid validate flag: " + arg);
            }
        }
        for(; argi < argc; argi++) {
            const char *path = argv[argi];

            PnaReader reader(path);
            for(uint64_t i  = 0; i < reader.getSequenceCount(); i++) {
                shared_ptr<PnaSequenceReader> seq = reader.openSequence(i, flags);
                printf(">%s %s\n",
                       seq->getMetadata().value("fasta.name"),
                       seq->getMetadata().value("fasta.comment"));
                       
                char buf[81];
                
                uint64_t rc;
                while( 0 != (rc = seq->read(buf, sizeof(buf) - 1)) ) {
                    buf[rc] = 0;
                    printf("%s\n", buf);
                }
            }
        }
    } else if((mode == "read") || (mode == "r")){
        bool packed = false;
        for(; argi < argc; argi++) {
            string arg = argv[argi];
            if(arg[0] != '-')
                break;

            if(arg == "--packed") {
                packed = true;
            } else {
                usage("Invalid validate flag: " + arg);
            }
        }
        for(; argi < argc; argi++) {
            const char *path = argv[argi];

            PnaReader reader(path);
            if(!packed) {
                for(uint64_t i  = 0; i < reader.getSequenceCount(); i++) {
                    shared_ptr<PnaSequenceReader> seq = reader.openSequence(i);
                    char *buf = new char[seq->size()];
                
                    seq->read(buf, seq->size());

                    delete [] buf;
                }
            } else {
                uint64_t buflen = reader.getMaxPackedBasesLength();
                uint8_t *buf = new uint8_t[buflen];
                PnaSequenceReader::packed_read_result_t result;

                for(uint64_t i  = 0; i < reader.getSequenceCount(); i++) {
                    shared_ptr<PnaSequenceReader> seq = reader.openSequence(i);
                    result = seq->packed_read(buf, buflen);
                }
            }
        }
    } else {
        usage("Invalid mode: '"+mode+"'");
    }
    
    return 0;
}

shared_ptr<PnaWriter> create_writer(const string &path) {
    shared_ptr<PnaWriter> fwriter(new PnaWriter(path.c_str()));
    create_file_metadata(fwriter);
    return fwriter;
}

void write_seq(FastaReader &reader,
               FastaSequenceDesc &seq,
               const string &path_fasta,
               int index_fasta,
               shared_ptr<PnaWriter> fwriter) {
    shared_ptr<PnaSequenceWriter> writer = fwriter->createSequence();

    char buf[1024 * 4];
    uint32_t rc;
    while( 0 != (rc = reader.read(buf, sizeof(buf))) ) {
        writer->write(buf, rc);
    }

    create_seq_metadata(writer, path_fasta, index_fasta, seq);
}

bool parse_ncbi_name(const string &name, map<string, string> &result) {
    vector<string> tokens = split(name, "|");
    if(tokens.empty() || (tokens.size() % 2 != 0))
        return false;

    for(size_t i = 0; i < tokens.size(); i += 2) {
        result[tokens[i]] = tokens[i+1];
    }

    return true;
}

bool parse_ncbi_comment(const string &comment, map<string, string> &result) {
    vector<string> tokens = split(comment, ",");
    if(tokens.size() < 2) {
        return false;
    }

    // species and molecule name
    vector<string> sm_tokens = split(tokens[0]);
    if(sm_tokens.size() < 3) {
        return false;
    }
    string species = join(sm_tokens.begin(), sm_tokens.begin() + 2);
    string molecule = join(sm_tokens.begin() + 2, sm_tokens.end());
    
    vector<string> assembly_tokens = split(tokens[1]);
    if(sm_tokens.size() < 1) {
        return false;
    }
    string assembly = join(assembly_tokens.begin(), assembly_tokens.end());

    result["species"] = species;
    result["molecule"] = molecule;
    result["assembly"] = assembly;

    return true;
}

void create_file_metadata(shared_ptr<PnaWriter> fwriter) {
    time_t rawtime;
    time(&rawtime);
    char *str = ctime(&rawtime);
    // strip trailing \n
    str[strlen(str) - 1] = 0;
    fwriter->addMetadata("create.time", str);
    str = getenv("USER");
    if(str)
        fwriter->addMetadata("create.user", str);
    str = getenv("HOSTNAME");
    if(str)
        fwriter->addMetadata("create.host", str);
}

void create_seq_metadata(shared_ptr<PnaSequenceWriter> writer,
                         const string &path_fasta,
                         int index_fasta,
                         FastaSequenceDesc &seq) {
    char strbuf[1024];

    sprintf(strbuf, "%d", index_fasta);
    writer->addMetadata("origin.index", strbuf);
    writer->addMetadata("origin.path", path_fasta.c_str());
    writer->addMetadata("origin.format", "fasta");
    writer->addMetadata("fasta.name", seq.name.c_str());
    writer->addMetadata("fasta.comment", seq.comment.c_str());

    {
        map<string,string> attrs;
        if(parse_ncbi_name(seq.name, attrs)) {
            string prefix = "ncbi.";
            for(auto kv: attrs) {
                writer->addMetadata((prefix+kv.first).c_str(),
                                    kv.second.c_str());
            }
        }
    }

    {
        map<string,string> attrs;
        if(parse_ncbi_comment(seq.comment, attrs)) {
            string prefix = "ncbi.";
            for(auto kv: attrs) {
                writer->addMetadata((prefix+kv.first).c_str(),
                                    kv.second.c_str());
            }
        }
    }

    {
        sprintf(strbuf, "%lu", writer->getBasePairCount());
        writer->addMetadata("size.base_pairs", strbuf);

        sprintf(strbuf, "%lu", writer->getByteCount());
        writer->addMetadata("size.bytes", strbuf);

        sprintf(strbuf, "%.2f%%",
                float(writer->getByteCount()) / writer->getBasePairCount() * 100);
        writer->addMetadata("size.compression", strbuf);
    }
}
