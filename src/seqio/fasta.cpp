#include "fasta.hpp"

#include "util.h"

#include <string.h>

using std::function;
using std::string;
using std::vector;
using namespace seqio::impl;

#define NCOLUMNS uint64_t(80)

namespace seqio {
    namespace impl {

        static vector<string> fasta_ext = {
            ".fasta", ".fas", ".fna", ".ffn", ".faa", ".frn", ".fa", ".mfa"
        };
        static vector<string> fasta_gzip_ext = {
            ".fasta.gz", ".fas.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".frn.gz", ".fa.gz", ".mfa.gz"
        };


        bool is_fasta_file_content(char const *path) {
            bool firstCol = true;
            char buf[4*1024];
            gzFile f;
            int rc;

            f = gzopen(path, "r");
            if(!f)
                goto fail;

            rc = gzread(f, buf, sizeof(buf));
            if(rc < 1)
                goto fail;

            for(int i = 0; i < rc; i++) {
                char c = buf[i];
                if(c == '>' && firstCol)
                    goto succeed;
                else if(c == '\n')
                    firstCol = true;
                else if(c <= 0)
                    goto fail;
                else
                    firstCol = false;
            }
            goto fail;

        succeed:
            gzclose(f);
            return true;

        fail:
            if(f) gzclose(f);
            return false;
        }

        bool is_fasta_file_name(char const *path) {
            for(string ext: fasta_ext)
                if(has_suffix(path, ext))
                    return true;

            return false;
        }

        bool is_fasta_gzip_file_name(char const *path) {
            for(string ext: fasta_gzip_ext)
                if(has_suffix(path, ext))
                    return true;

            return false;
        }

    }
}
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
                switch(base) {
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

bool FastaMetadata::hasKey(char const *key) const {
    return (0 == strcmp(key, SEQIO_KEY_NAME))
        || (0 == strcmp(key, SEQIO_KEY_COMMENT));
            
}

uint32_t FastaMetadata::getKeyCount() const {
    return 2;
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

IConstDictionary const &FastaSequence::getMetadata() {
    return metadata;
}

uint64_t FastaSequence::read(char *buffer,
                             uint64_t buffer_length) {
    if(parse.eos)
        return 0;

    uint64_t n = 0;
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
 * CLASS FastaWriter
 *
 **********************************************************************/

FastaWriter::FastaWriter(char const *path,
                         seqio_file_format file_format) {
    inSequence = false;

    switch(file_format) {
    case SEQIO_FILE_FORMAT_FASTA: {
        FILE *f = fopen(path, "w");

        doWrite = [=] (char const *buffer, uint32_t len) {
            size_t rc = fwrite(buffer, 1, len, f);
            if(rc != len)
                raise_io("Failed writing to %s", path);
        };

        doClose = [=] () {
            fclose(f);
        };
    } break;
    case SEQIO_FILE_FORMAT_FASTA_GZIP: {
        gzFile f = gzopen(path, "w");

        doWrite = [=] (char const *buffer, uint32_t len) {
            size_t rc = gzwrite(f, buffer, len);
            if(rc != len)
                raise_io("Failed writing to %s", path);
        };

        doClose = [=] () {
            gzclose(f);
        };
    } break;
    default:
        panic();
    }
}

FastaWriter::~FastaWriter() {
    doClose();
}

void FastaWriter::createSequence(IConstDictionary const *metadata) {
    name = metadata->getValue(SEQIO_KEY_NAME);
    if(metadata->hasKey(SEQIO_KEY_COMMENT))
        comment = metadata->getValue(SEQIO_KEY_COMMENT);

    inSequence = true;
    writeMetadata();
    column = 0;
}

void FastaWriter::write(char const *buffer,
                        uint64_t length) {
    if(!inSequence) {
        raise_state("Must create sequence");
    }

    while(length > 0) {
        uint64_t write_length = std::min(NCOLUMNS - column, length);
        doWrite(buffer, write_length);
        doWrite("\n", 1);
        
        buffer += write_length;
        length -= write_length;

        column += write_length;
        if(column == NCOLUMNS)
            column = 0;
    }
}

void FastaWriter::writeMetadata() {
    if(name.length() == 0)
        raise_state("Must specify sequence name");

    doWrite(">", 1);
    doWrite(name.c_str(), name.length());

    if(comment.length() > 0) {
        doWrite(" ", 1);
        doWrite(comment.c_str(), comment.length());
    }
    
    doWrite("\n", 1);
}
