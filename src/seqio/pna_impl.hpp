#pragma once

#include "pna.hpp"
#include "seqio_impl.hpp"

namespace seqio {
    namespace impl {

        bool is_pna_file_content(char const *path);
        bool is_pna_file_name(char const *path);

/**********************************************************************
 *
 * CLASS PnaMetadata
 *
 **********************************************************************/
        class PnaMetadata : public IConstDictionary {
        public:
            PnaMetadata(pna::PnaMetadata const &metadata);
            virtual ~PnaMetadata();

            virtual uint32_t getKeyCount() const override;
            virtual char const *getKey(uint32_t key_index) const override;
            virtual char const *getValue(char const *key) const override;

        private:
            pna::PnaMetadata const metadata;
        };

/**********************************************************************
 *
 * CLASS PnaSequence
 *
 **********************************************************************/
        class PnaSequence : public ISequence {
        public:
            PnaSequence(std::shared_ptr<pna::PnaSequenceReader> reader_);
            virtual ~PnaSequence();

            virtual IConstDictionary const &getMetadata() override;
            virtual uint32_t read(char *buffer,
                                  uint32_t buffer_length) override;

        private:
            std::shared_ptr<pna::PnaSequenceReader> reader;
            PnaMetadata metadata;
        };

/**********************************************************************
 *
 * CLASS PnaSequenceIterator
 *
 **********************************************************************/
        class PnaSequenceIterator : public ISequenceIterator {
        public:
            PnaSequenceIterator(char const *path,
                                seqio_base_transform transform);
            virtual ~PnaSequenceIterator();

            virtual ISequence *nextSequence() override;

        private:
            std::shared_ptr<pna::PnaReader> reader;
            uint64_t index;
        };

/**********************************************************************
 *
 * CLASS PnaWriter
 *
 **********************************************************************/
        class PnaWriter : public IWriter {
        public:
            PnaWriter(char const *path,
                      seqio_file_format file_format);
            virtual ~PnaWriter();

            virtual void createSequence() override;
            virtual void addMetadata(char const *key,
                                     char const *value) override;
            virtual void write(char const *buffer,
                               uint32_t length) override;

        private:
            std::shared_ptr<pna::PnaWriter> writer;
            std::shared_ptr<pna::PnaSequenceWriter> sequence;
        };

    }
}
