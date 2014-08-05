#include "fasta.hpp"

#include "util.h"

#include <string.h>

using std::function;
using std::string;
using std::vector;
using namespace seqio::impl;

/**********************************************************************
 *
 * CLASS FileHandle
 *
 **********************************************************************/
__FileHandle::__FileHandle(gzFile f_)
    : f(f_) {
}

__FileHandle::~__FileHandle() {
    gzclose(f);
    f = 0;
}

__FileHandle::operator gzFile() {
    return f;
}

/**********************************************************************
 *
 * CLASS FastaRawStream
 *
 **********************************************************************/
FastaRawStream::FastaRawStream(FileHandle f_,
                               z_off_t start_)
    : f(f_)
    , start(start_) {

    fstate.read_offset = start;
    cache.len = 0;
    cache.index = 0;
}

FastaRawStream *FastaRawStream::createSubstream() {
    FastaRawStream *substream = new FastaRawStream;
    substream->cache = cache;
    substream->fstate = fstate;
    substream->f = f;
    substream->start = tell_abs();
    return substream;
}

int FastaRawStream::nextChar() {
    if(fstate.eof) {
        return -1;
    }

    if(cache.index == cache.len) {
        z_off_t read_offset = fstate.read_offset + cache.len;
        z_off_t off = gzseek(*f, read_offset, SEEK_SET);
        if(off != read_offset) raise_io("Failed seeking to offset %zu", (size_t)off);

        int rc = gzread(*f, cache.buf, sizeof(cache.buf));
        if(rc < 0) {
            raise(IO, "Failed reading from file.");
        }
        if(rc == 0) {
            fstate.eof = true;
            return -1;
        }
        cache.index = 0;
        cache.len = rc;
        fstate.read_offset = read_offset;
    }

    return (int)cache.buf[cache.index++];
}

z_off_t FastaRawStream::tell_abs() {
    return fstate.read_offset + cache.index;
}

void FastaRawStream::seek_abs(z_off_t offset) {
    if(offset != gzseek(*f, offset, SEEK_SET)) raise_io("Failed seeking");
    fstate.read_offset  = offset;
    cache.len = 0;
    cache.index = 0;
}

/**********************************************************************
 *
 * CLASS CharInterpreter
 *
 **********************************************************************/
CharInterpreter::CharInterpreter(seqio_base_transform transform) {
    for(int c = 0; c < 256; c++) {
        char base = c;
        Action firstCol;
        Action otherCol;

        if(c == '\n') {
            firstCol = otherCol = NEWLINE;
        } else if(!isgraph(c)) {
            firstCol = otherCol = IGNORE;
        } else if(c == '>') {
            firstCol = HEADER;
            otherCol = APPEND_SEQUENCE;
        } else {
            firstCol = otherCol = APPEND_SEQUENCE;
            if(transform == SEQIO_BASE_TRANSFORM_CAPS_GATCN) {
                base = toupper(c);
                switch(c) {
                case 'G':
                case 'A':
                case 'T':
                case 'C':
                    // no-op
                    break;
                default:
                    base = 'N';
                    break;
                }
            }
        }

        bases[c] = base;
        actions_firstCol[c] = firstCol;
        actions_otherCol[c] = otherCol;
    }
}

/**********************************************************************
 *
 * CLASS FastaMetadata
 *
 **********************************************************************/
FastaMetadata::FastaMetadata(string name_, string comment_)
    : name(name_)
    , comment(comment_) {
}

FastaMetadata::~FastaMetadata() {
}

uint32_t FastaMetadata::getKeyCount() const {
    return 1;
}

char const *FastaMetadata::getKey(uint32_t key_index) const {
    switch(key_index) {
    case 0:
        return SEQIO_KEY_NAME;
    case 1:
        return SEQIO_KEY_COMMENT;
    default:
        raise_parm("Requested key index %zu exceeds key count 2", size_t(key_index));
    }
}

char const *FastaMetadata::getValue(char const *key) const {
    if(0 == strcmp(key, SEQIO_KEY_NAME)) {
        return name.c_str();
    } else if(0 == strcmp(key, SEQIO_KEY_COMMENT)) {
        return comment.c_str();
    } else {
        raise_parm("Unknown key: %s", key);
    }
}

/**********************************************************************
 *
 * CLASS FastaSequence
 *
 **********************************************************************/
FastaSequence::FastaSequence(FastaMetadata const &metadata_,
                             FastaRawStream *stream_,
                             CharInterpreter const *interpreter_,
                             function<void (FastaSequence *sequence)> onClose_)
    : metadata(metadata_)
    , stream(stream_)
    , interpreter(interpreter_)
    , onClose(onClose_) {

    parse.firstCol = true;
    parse.eos = false;
}

FastaSequence::~FastaSequence() {
    onClose(this);
    delete stream;
}

IMetadata const &FastaSequence::getMetadata() {
    return metadata;
}

uint32_t FastaSequence::read(char *buffer,
                             uint32_t buffer_length) {
    if(parse.eos)
        return 0;

    uint32_t n = 0;
    int c;

    while((n < buffer_length) && (c = stream->nextChar()) != -1) {
        CharInterpreter::Action action = interpreter->getAction(c, parse.firstCol);

        switch(action) {
        case CharInterpreter::IGNORE:
            parse.firstCol = false;
            break;
        case CharInterpreter::NEWLINE:
            parse.firstCol = true;
            break;
        case CharInterpreter::APPEND_SEQUENCE:
            parse.firstCol = false;
            buffer[n++] = interpreter->getBase(c);
            break;
        case CharInterpreter::HEADER:
            parse.eos = true;
            parse.eos_offset = stream->tell_abs() - 1;
            goto break_outer;
        default:
            panic();
        }
    }
break_outer:

    if(c == -1) {
        parse.eos = true;
        parse.eos_offset = stream->tell_abs();
    }

    return n;
}

z_off_t FastaSequence::tellEnd() {
    if(!parse.eos) {
        z_off_t offset = stream->tell_abs();

        char buf[1024];
        while(!parse.eos)
            read(buf, sizeof(buf));

        parse.eos = false;
        stream->seek_abs(offset);
    }

    return parse.eos_offset;
}

/**********************************************************************
 *
 * CLASS FastaSequenceIterator
 *
 **********************************************************************/
FastaSequenceIterator::FastaSequenceIterator(char const *path,
                                             seqio_base_transform transform)
    : interpreter(transform) {

    callback = std::make_shared<Callback>(this);

    FileHandle f = std::make_shared<__FileHandle>( gzopen(path, "r") );
    if(!f->f) raise_io("Failed opening %s", path);
    
    stream = new FastaRawStream(f, 0);
    currSequence = nullptr;
    eos_offset = 0;
}

FastaSequenceIterator::~FastaSequenceIterator() {
    callback->iteratorClosing();

    delete stream;
}


ISequence *FastaSequenceIterator::nextSequence() {
    if(currSequence) {
        eos_offset = currSequence->tellEnd();
        currSequence = nullptr;
    }

    stream->seek_abs(eos_offset);

    // ---
    // --- Find next header
    // ---
    {
        bool firstCol = true;
        bool foundHeader = false;
        int c;
        while((c = stream->nextChar()) > -1) {
            if(c == '\n') {
                firstCol = true;
            } else if((c == '>') && firstCol) {
                foundHeader = true;
                break;
            } else {
                firstCol = false;
            }
        }

        if(!foundHeader) return nullptr;
    }

    // ---
    // --- Get name and comment
    // ---
    string name;
    string comment;
    {
        int c;
        while(!isspace(c = stream->nextChar()) && (c > -1)) {
            name += c;
        }
        if(c == -1) return nullptr;

        if((c != '\n') && (c != '\r')) {
            while(((c = stream->nextChar()) != '\n') && (c > -1)) {
                if(c != '\r')
                    comment += c;
            }
            if(c == -1) return nullptr;
        }
    }

    // The functor will maintain a shared pointer to the callback, meaning
    // the callback will still be valid even if the iterator has been disposed.
    std::shared_ptr<Callback> callback_ = callback;
    auto onClose = [callback_] (FastaSequence *sequence) {
        callback_->sequenceClosing(sequence);
    };

    currSequence = new FastaSequence(FastaMetadata(name, comment),
                                     stream->createSubstream(),
                                     &interpreter,
                                     onClose);

    return currSequence;
}

FastaSequenceIterator::Callback::Callback(FastaSequenceIterator *thiz_)
    : thiz(thiz_) {
}

void FastaSequenceIterator::Callback::sequenceClosing(FastaSequence *sequence) {
    if(thiz && (sequence == thiz->currSequence)) {
        thiz->eos_offset = sequence->tellEnd();
        thiz->currSequence = nullptr;
    }
}

void FastaSequenceIterator::Callback::iteratorClosing() {
    thiz = nullptr;
}

/**********************************************************************
 *
 * OBSOLETE
 *
 **********************************************************************/
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
    abort();
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
