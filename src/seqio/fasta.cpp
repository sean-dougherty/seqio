#include "fasta.h"

#include "seqio_impl.hpp"
#include "util.h"

#include <string.h>

using namespace seqio;

FastaRawStream::FastaRawStream(gzFile f_, seqio_base_transform transform)
    : f(f_) {

    fstate.offset = gztell(f);
    start = fstate.offset;

    initActions(transform);
}

FastaRawStream *FastaRawStream::createSubstream() {
    FastaRawStream *substream = new FastaRawStream;
    memcpy(substream, this, sizeof(*this));
    substream->start = tell();
    return substream;
}

int FastaRawStream::nextChar() {
    z_off_t off = gzseek(f, fstate.offset, SEEK_SET);
    if(off != fstate.offset) raise(IO, "Failed seeking to offset %zu", (size_t)off);

    int rc = gzread(f, cache.buf, sizeof(cache.buf));
    if(rc < 0) {
        raise(IO, "Failed reading from file.");
    } else if(rc == 0) {
        fstate.eof = true;
        return -1;
    }
    cache.index = 0;
    cache.len = rc;
    fstate.offset += rc;

    return (int)cache.buf[cache.index++];
}

z_off_t FastaRawStream::tell() {
    return fstate.offset + cache.index;
}

void FastaRawStream::initActions(seqio_base_transform transform) {
    for(int c = 0; c < 256; c++) {
        CharActionEntry &entry = actions[c];

        entry.c = c;

        if((c == '\r') || (c == '\n')) {
            entry.firstCol = entry.otherCol = CharActionEntry::NEWLINE;
        } else if(!isgraph(c)) {
            entry.firstCol = entry.otherCol = CharActionEntry::IGNORE;
        } else if(c == '>') {
            entry.firstCol = CharActionEntry::HEADER;
            entry.otherCol = CharActionEntry::APPEND;
        } else {
            entry.firstCol = entry.otherCol = CharActionEntry::APPEND;
            if(transform == SEQIO_BASE_TRANSFORM_CAPS_GATCN) {
                entry.c = toupper(entry.c);
                switch(entry.c) {
                case 'G':
                case 'A':
                case 'T':
                case 'C':
                    // no-op
                    break;
                default:
                    entry.c = 'N';
                    break;
                }
            }
        }
    }
}

FastaReader::FastaReader(const char *path, Translate translate) {
    errif(NULL == (f = gzopen(path, "r")),
          "Failed opening %s", path);

    initActions(translate);
}

FastaReader::~FastaReader() {
    close();
}

bool FastaReader::nextSequence(FastaSequenceDesc &desc) {
    errif(state == SEQ, "Currently in sequence");

    desc.name.clear();
    desc.comment.clear();

    if(state == END) return false;

    int c;
    if(state == INIT) {
        bool foundHeader = false;
        while((c = nextChar()) > -1) {
            if((c == '\n') || (c == '\r')) {
                firstCol = true;
            } else if((c == '>') && firstCol) {
                foundHeader = true;
                break;
            } else {
                firstCol = false;
            }
        }

        if(!foundHeader) return false;
    }

    while(!isspace(c = nextChar()) && (c > -1)) {
        desc.name += c;
    }

    if(c == -1) return false;

    if((c != '\n') && (c != '\r')) {
        while(((c = nextChar()) != '\n') && (c != '\r') && (c > -1)) {
            desc.comment += c;
        }
        if(c == -1) return false;
    }

    state = SEQ;
    firstCol = true;

    return true;
}

uint32_t FastaReader::read(char *buf, uint32_t buflen) {
    if(state != SEQ) return 0;

    uint32_t n = 0;
    int c;

    while((n < buflen) && (c = nextChar()) != -1) {
        CharActionEntry entry = actions[c];

        if(firstCol) {
            CharActionEntry::Action action = entry.firstCol;
            switch(action) {
            case CharActionEntry::APPEND:
                firstCol = false;
                buf[n++] = entry.c;
                break;
            case CharActionEntry::NEWLINE:
                // no-op
                break;
            case CharActionEntry::HEADER:
                firstCol = false;
                state = HEADER;
                goto break_outer;
                break;
            default:
                err("Program logic error");
                break;
            }
        } else {
            CharActionEntry::Action action = entry.otherCol;
            switch(action) {
            case CharActionEntry::APPEND:
                buf[n++] = entry.c;
                break;
            case CharActionEntry::NEWLINE:
                firstCol = true;
                break;
            default:
                err("Program logic error");
            }
        }            
    }
break_outer:

    return n;
}

void FastaReader::close() {
    if(!f) return;

    errif(0 > gzclose(f),
          "Failed closing gzip file");
    f = NULL;
}

int FastaReader::nextChar() {
    throw 1;
}

void FastaReader::initActions(Translate translate) {
    for(int c = 0; c < 256; c++) {
        CharActionEntry &entry = actions[c];

        entry.c = c;

        if((c == '\r') || (c == '\n')) {
            entry.firstCol = entry.otherCol = CharActionEntry::NEWLINE;
        } else if(!isgraph(c)) {
            entry.firstCol = entry.otherCol = CharActionEntry::IGNORE;
        } else if(c == '>') {
            entry.firstCol = CharActionEntry::HEADER;
            entry.otherCol = CharActionEntry::APPEND;
        } else {
            entry.firstCol = entry.otherCol = CharActionEntry::APPEND;
            if(translate == Translate_Caps_GATCN) {
                entry.c = toupper(entry.c);
                switch(entry.c) {
                case 'G':
                case 'A':
                case 'T':
                case 'C':
                    // no-op
                    break;
                default:
                    entry.c = 'N';
                    break;
                }
            }
        }
    }
}
