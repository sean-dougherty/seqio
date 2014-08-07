#pragma once

#include "seqio_impl.hpp"

#include <ctype.h>
#include <stdint.h>
#include <zlib.h>

#include <memory>
#include <string>
#include <vector>

namespace seqio {
    namespace impl {

        bool is_fasta_file_content(char const *path);
        bool is_fasta_file_name(char const *path);
        bool is_fasta_gzip_file_name(char const *path);

/**********************************************************************
 *
 * CLASS FileHandle
 *
 **********************************************************************/
        struct __FileHandle {
            __FileHandle(gzFile f_);
            ~__FileHandle();

            operator gzFile();

            gzFile f;
        };

        typedef std::shared_ptr<__FileHandle> FileHandle;

/**********************************************************************
 *
 * CLASS FastaRawStream
 *
 **********************************************************************/
        class FastaRawStream {
        public:
            FastaRawStream(FileHandle f_,
                           z_off_t start_);

            int nextChar();
            z_off_t tell_abs();
            void seek_abs(z_off_t offset);

            FastaRawStream *createSubstream();

        private:
            FastaRawStream() {}

            struct {
                char buf[1024*64];
                uint32_t len;
                uint32_t index;
            } cache;

            struct {
                z_off_t read_offset;
                bool eof = false;
            } fstate;

            FileHandle f;
            z_off_t start;
        };

/**********************************************************************
 *
 * CLASS CharInterpreter
 *
 **********************************************************************/
        class CharInterpreter {
        public:
            enum Action {
                IGNORE,
                NEWLINE,
                APPEND_SEQUENCE,
                HEADER
            };

            CharInterpreter(seqio_base_transform transform);

            inline char getBase(char c) const {
                return bases[c & 0xFF];
            }

            inline Action getAction(char c, bool firstCol) const {
                if(firstCol) 
                    return actions_firstCol[c & 0xFF];
                else
                    return actions_otherCol[c & 0xFF];
            }

        private:
            char bases[256];
            Action actions_firstCol[256];
            Action actions_otherCol[256];
        };

/**********************************************************************
 *
 * CLASS FastaMetadata
 *
 **********************************************************************/
        class FastaMetadata : public IConstDictionary {
        public:
            FastaMetadata(std::string name, std::string comment);
            virtual ~FastaMetadata();

            virtual bool hasKey(char const *key) const override;
            virtual uint32_t getKeyCount() const override;
            virtual char const *getKey(uint32_t key_index) const override;
            virtual char const *getValue(char const *key) const override;

        private:
            std::string name;
            std::string comment;
        };

/**********************************************************************
 *
 * CLASS FastaSequence
 *
 **********************************************************************/
        class FastaSequence : public ISequence {
        public:
            FastaSequence(FastaMetadata const &metadata_,
                          FastaRawStream *stream_,
                          CharInterpreter const *interpreter_,
                          std::function<void (FastaSequence *sequence)> onClose_);
            virtual ~FastaSequence();

            virtual IConstDictionary const &getMetadata() override;
            virtual uint32_t read(char *buffer,
                                  uint32_t buffer_length) override;
            z_off_t tellEnd();

        private:
            FastaMetadata metadata;
            FastaRawStream * const stream;
            CharInterpreter const * const interpreter;
            std::function<void (FastaSequence *sequence)> onClose;
            struct {
                bool firstCol;
                bool eos;
                z_off_t eos_offset;
            } parse;
        };

/**********************************************************************
 *
 * CLASS FastaSequenceIterator
 *
 **********************************************************************/
        class FastaSequenceIterator : public ISequenceIterator {
        public:
            FastaSequenceIterator(char const *path,
                                  seqio_base_transform transform);
            virtual ~FastaSequenceIterator();

            virtual ISequence *nextSequence() override;

        private:
            class Callback {
            public:
                Callback(FastaSequenceIterator *thiz_);

                void sequenceClosing(FastaSequence *sequence);
                void iteratorClosing();

            private:
                FastaSequenceIterator *thiz;
            };
            std::shared_ptr<Callback> callback;

            FastaRawStream *stream;
            CharInterpreter interpreter;
            FastaSequence *currSequence;
            z_off_t eos_offset;
        };

/**********************************************************************
 *
 * CLASS FastaWriter
 *
 **********************************************************************/
        class FastaWriter : public IWriter {
        public:
            FastaWriter(char const *path,
                        seqio_file_format file_format);
            virtual ~FastaWriter();

            virtual void createSequence(IConstDictionary const *metadata) override;
            virtual void write(char const *buffer,
                               uint32_t length) override;

        private:
            void writeMetadata();

            std::function<void (char const *buffer, uint32_t length)> doWrite;
            std::function<void ()> doClose;

            bool inSequence;
            std::string name;
            std::string comment;
            int column;
        };
    }
}
