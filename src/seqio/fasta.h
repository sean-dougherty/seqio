#pragma once

#include "seqio.h"

#include <ctype.h>
#include <stdint.h>
#include <zlib.h>

#include <string>

namespace seqio {

    struct FastaSequenceDesc {
        std::string name;
        std::string comment;
    };

    class FastaRawStream {
    public:
        FastaRawStream(gzFile file, seqio_base_transform transform);

        int nextChar();
        z_off_t tell();
        void reset();

        FastaRawStream *createSubstream();

    private:
        FastaRawStream() {}

        void initActions(seqio_base_transform transform);

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

        struct {
            char buf[1024*16];
            uint16_t len = 0;
            uint16_t index = 0;
        } cache;

        struct {
            z_off_t offset;
            bool eof = false;
        } fstate;

        gzFile f;
        z_off_t start;
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
