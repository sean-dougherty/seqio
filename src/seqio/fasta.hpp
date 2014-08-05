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
                uint16_t len;
                uint16_t index;
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
        class FastaMetadata : public IMetadata {
        public:
            FastaMetadata(std::string name, std::string comment);
            virtual ~FastaMetadata();

            virtual uint32_t getKeyCount() const;
            virtual char const *getKey(uint32_t key_index) const;
            virtual char const *getValue(char const *key) const;

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

            virtual IMetadata const &getMetadata();
            virtual uint32_t read(char *buffer,
                                  uint32_t buffer_length);
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

            virtual ISequence *nextSequence();

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
 * OBSOLETE
 *
 **********************************************************************/
        struct FastaSequenceDesc {
            std::string name;
            std::string comment;
        };

        class FastaReader {
        public:
            enum Translate {
                Translate_None,
                Translate_Caps_GATCN
            };

            FastaReader(const char *path, Translate translate = Translate_None);
            ~FastaReader();

            bool nextSequence(FastaSequenceDesc &desc);
            uint32_t read(char *buf, uint32_t buflen);
            void close();

        private:
            struct CharActionEntry {
                enum Action {
                    IGNORE,
                    NEWLINE,
                    APPEND,
                    HEADER
                };

                Action firstCol;
                Action otherCol;
                char c;
            } actions[256];

            void initActions(Translate translate);

            int nextChar();
            enum {
                INIT, SEQ, HEADER, END
            } state = INIT;
            bool firstCol = true;
            struct {
                char buf[1024*4];
                uint16_t len = 0;
                uint16_t index = 0;
            } cache;
            gzFile f;
        };

    }
}
