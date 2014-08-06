#pragma once

#include "pna_layout.h"

#include <stdint.h>
#include <stdio.h>

#include <atomic>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace seqio {
    namespace pna {

        bool is_pna_file_content(char const *path);
        bool is_pna_file_name(char const *path);

        // todo: confirm can be same
        const uint16_t READBUF_CAPACITY = 4 * 1024;
        const uint16_t WRITEBUF_CAPACITY = 16 * 1024;
    
        enum base_t {A = 0, C = 1, G = 2, T = 3, N = 4};

        class PnaMetadata {
            friend class PnaReader;

            PnaMetadata();
            PnaMetadata(const metadata_t &metadata,
                        uint8_t *file_start,
                        const char *strings);
        public:
            uint32_t size() const;
            bool pair(uint32_t index, const char **key, const char **value) const;
            const char *value(const char *key) const;

        private:
            struct {
                metadata_entry_t *begin;
                metadata_entry_t *end;
                uint32_t size;
            } entries;
            const char *strings;
        };

        class FilePointerPool;

        class FilePointerGuard {
            friend class FilePointerPool;

            FilePointerGuard(FilePointerPool *pool, FILE *f);
        public:
            ~FilePointerGuard();

            void manage(FILE **field);
            void manual_release();

        private:
            FilePointerPool *pool;
            FILE *f;
            FILE **f_managed = nullptr;
        };

        class FilePointerPool {
        public:
            FilePointerPool(const std::string &path);
            ~FilePointerPool();

            FilePointerGuard acquire();

        private:
            friend class FilePointerGuard;
            void release(FILE *);

            std::string path;
            std::atomic_flag lock = ATOMIC_FLAG_INIT;
            std::vector<FILE *> handles;
        };

        class PnaSequenceReader {
        public:
            enum Flags {
                Standard = 0,
                IgnoreN = (1 << 0)
            };
        private:
            friend class PnaReader;

            PnaSequenceReader(FilePointerGuard fguard,
                              const sequence_t &sequence,
                              const PnaMetadata &metadata,
                              uint32_t flags);

        public:
            ~PnaSequenceReader();

            uint64_t size();
            void seek(uint64_t offset);
            uint64_t read(char *buf, uint64_t buflen);
            const PnaMetadata getMetadata();
            void close();

            struct packed_read_result_t {
                uint64_t bases_count;
                struct {
                    const seqfragment_t *begin;
                    uint64_t count;
                } seqfragments;
                struct {
                    uint8_t *begin;
                    uint64_t buflen;
                    uint64_t count;
                } packed_bases;
            };

            packed_read_result_t packed_read(uint8_t *packed_buf, uint64_t buflen);

        private:
            seqfragment_t *find_next_seqfragment(uint64_t offset);

            FilePointerGuard fguard;
            FILE *fpna;
            sequence_t sequence;
            PnaMetadata metadata;
            uint32_t flags;
            struct {
                seqfragment_t *begin;
                seqfragment_t *next;
                seqfragment_t *end;
            } seqfragments;
            uint8_t shift = 0;
            uint64_t seqOffset = 0;
            uint32_t *packedByteLookup;
            struct {
                unsigned char *buf;
                uint16_t len = 0;
                uint16_t index = 0;
                uint64_t bases_offset = 0;
                unsigned char curr;
            } packedCache;
        };

        class PnaReader {
        public:
            PnaReader(const char *path);
            ~PnaReader();

            uint64_t getSequenceCount();
            uint64_t getMaxSequenceFragments();
            uint64_t getMaxPackedBasesLength();
            const PnaMetadata getMetadata();
            const PnaMetadata getSequenceMetadata(uint64_t index);
            std::shared_ptr<PnaSequenceReader> openSequence(uint64_t index,
                                                            uint32_t flags = PnaSequenceReader::Standard);

        private:
            std::string path;
            FilePointerPool fpool;
            header_t header;
            const sequence_t *sequences;
            PnaMetadata metadata;
            const char *strings;
            struct {
                void *addr = nullptr;
                uint64_t length;
                uint8_t *file_start;
            } mmap;
        };

// todo: move inside class
        typedef uint32_t stringid_t;

        class StringStorageWriter {
        public:
            stringid_t getId(const char *str);
            void write(FILE *f, string_storage_t &header);
            uint32_t getOffset(stringid_t id);
    
        private:
            typedef std::map<std::string, stringid_t> String2Id;
            String2Id ids;
            typedef std::map<stringid_t, uint32_t> OffsetMap;
            OffsetMap offsets;
        };

// todo: move inside writer
        typedef std::map<stringid_t, stringid_t> StringIdMap;

        class MetadataWriter {
        public:
            MetadataWriter(metadata_t &metadata,
                           StringStorageWriter &strings);

            void addMetadata(const char *key, const char *value);
            void write(FILE *f);

        private:
            metadata_t &metadata;
            StringStorageWriter &strings;
            StringIdMap idmap;
        };

        class PnaSequenceWriter {
            friend class PnaWriter;

            PnaSequenceWriter(FILE *f,
                              sequence_t &sequence,
                              MetadataWriter &metadata);

        public:
            ~PnaSequenceWriter();

            void write(char const *buf, uint64_t buflen);
            void addMetadata(const char *key, const char *value);
            uint64_t getBasePairCount();
            uint64_t getByteCount();
            void close();

        private:
            void addPackedByte();
            void flushCache();

            FILE *fpna;
            base_t base_map[256];
            uint64_t seqOffset = 0;
            unsigned char packedByte = 0;
            int shift = 0;
            bool in_seqfragment = false;
            seqfragment_t seqfragment;
            std::vector<seqfragment_t> seqfragments;
            sequence_t &sequence;
            struct {
                char buf[WRITEBUF_CAPACITY];
                uint16_t len = 0;
            } packedCache;
            MetadataWriter &metadata;
        };

        class PnaWriter {
        public:
            PnaWriter(const char *path);
            ~PnaWriter();

            void addMetadata(const char *key, const char *value);

            std::shared_ptr<PnaSequenceWriter> createSequence();

            void close();

        private:
            void closeSequence();

            FILE *f = nullptr;
            header_t header;
            std::vector<sequence_t *> sequences;
            StringStorageWriter strings;
            MetadataWriter metadata;
            std::vector<MetadataWriter *> sequences_metadata;
            std::weak_ptr<PnaSequenceWriter> activeSequenceWriter;
        };
    }
}
