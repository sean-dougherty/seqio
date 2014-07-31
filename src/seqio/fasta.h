#pragma once

#include <ctype.h>
#include <stdint.h>
#include <zlib.h>

#include <string>

namespace seqio {

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
