#include "pna_impl.hpp"

using namespace seqio;
using namespace seqio::impl;

namespace seqio {
    namespace impl {

        bool is_pna_file_content(char const *path) {
            return pna::is_pna_file_content(path);
        }

        bool is_pna_file_name(char const *path) {
            return pna::is_pna_file_name(path);
        }

    }
}

/**********************************************************************
 *
 * CLASS PnaMetadata
 *
 **********************************************************************/
PnaMetadata::PnaMetadata(pna::PnaMetadata const &metadata_)
    : metadata(metadata_) {
}

PnaMetadata::~PnaMetadata() {
}

bool PnaMetadata::hasKey(char const *key) const {
    return metadata.value(key) != nullptr;
}

uint32_t PnaMetadata::getKeyCount() const {
    return metadata.size();
}

char const *PnaMetadata::getKey(uint32_t key_index) const {
    char const *key, *value;

    if(!metadata.pair(key_index, &key, &value)) {
        raise_parm("Invalid key index: %zu", size_t(key_index));
    }

    return key;
}

char const *PnaMetadata::getValue(char const *key) const {
    char const *value =  metadata.value(key);

    if(!value) {
        raise_parm("Invalid key: %s", key);
    }

    return value;
}

/**********************************************************************
 *
 * CLASS PnaSequence
 *
 **********************************************************************/
PnaSequence::PnaSequence(std::shared_ptr<pna::PnaSequenceReader> reader_)
    : reader(reader_)
    , metadata(reader_->getMetadata()) {
}

PnaSequence::~PnaSequence() {
}

IConstDictionary const &PnaSequence::getMetadata() {
    return metadata;
}

uint64_t PnaSequence::read(char *buffer,
                           uint64_t buffer_length) {
    return reader->read(buffer, buffer_length);
}

/**********************************************************************
 *
 * CLASS PnaSequenceIterator
 *
 **********************************************************************/
PnaSequenceIterator::PnaSequenceIterator(char const *path,
                                         seqio_base_transform transform)
    : reader(std::make_shared<pna::PnaReader>(path))
    , index(0) {
}

PnaSequenceIterator::~PnaSequenceIterator() {
}

ISequence *PnaSequenceIterator::nextSequence() {
    if(index >= reader->getSequenceCount()) return nullptr;

    PnaSequence *sequence = new PnaSequence(reader->openSequence(index++));
    return sequence;
}

/**********************************************************************
 *
 * CLASS PnaWriter
 *
 **********************************************************************/
PnaWriter::PnaWriter(char const *path,
                     seqio_file_format file_format)
    : writer(std::make_shared<pna::PnaWriter>(path)) {
}

PnaWriter::~PnaWriter() {
}

void PnaWriter::createSequence(IConstDictionary const *metadata) {
    sequence = writer->createSequence();
    for(uint32_t i = 0; i < metadata->getKeyCount(); i++) {
        char const *key = metadata->getKey(i);
        char const *value = metadata->getValue(key);
        sequence->addMetadata(key, value);
    }
}

void PnaWriter::write(char const *buffer,
                      uint64_t length) {
    sequence->write(buffer, length);
}
