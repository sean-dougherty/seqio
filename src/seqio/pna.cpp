#define _FILE_OFFSET_BITS   64

//#define DEBUG_NBLOCK

#include <string.h>
#include <unistd.h>
#include <sys/mman.h>

#include <algorithm>

#include "pna.hpp"
#include "seqio_impl.hpp"
#include "util.h"

using namespace std;
using namespace seqio::pna;

#define MAX_SEQFRAGMENT_LEN ((uint32_t)~0)
#define MAX_STRING_STORAGE ((uint32_t)~0)

namespace seqio {
    namespace pna {

        bool is_pna_file_content(char const *path) {
            FILE *f = fopen(path, "r");

            if(!f)
                goto fail;
    
            header_t header;
            if(1 != fread(&header, sizeof(header), 1, f))
                goto fail;

            if(header.signature != PNA_FILE_SIGNATURE)
                goto fail;

            fclose(f);
            return true;

        fail:
            if(f) fclose(f);
            return false;
        }

        bool is_pna_file_name(char const *path) {
            return has_suffix(path, ".pna");
        }

    }
}

PnaMetadata::PnaMetadata() 
    : entries({nullptr, nullptr, 0l})
    , strings(nullptr)
{
}

PnaMetadata::PnaMetadata(const metadata_t &metadata,
                         uint8_t *file_start,
                         const char *strings_) {
    entries.begin = (metadata_entry_t *)(file_start + metadata.entries_filepos);
    entries.size = metadata.entries_count;
    entries.end = entries.begin  + entries.size;
    strings = strings_;
}

uint32_t PnaMetadata::size() const {
    return entries.size;
}

bool PnaMetadata::pair(uint32_t index,
                       const char **key,
                       const char **value) const {
    if(index >= size())
        return false;

    metadata_entry_t *entry = entries.begin + index;
    *key = strings + entry->key;
    *value = strings + entry->value;

    return true;
}

const char *PnaMetadata::value(const char *key) const {
    struct compare_t {
        const char *strings;
        const char *key;

        static int cmp(const void *a, const void *b) {
            const compare_t *key = (const compare_t *)a;
            const metadata_entry_t *entry = (const metadata_entry_t *)b;

            return strcmp(key->key, key->strings + entry->key);
        }
    } compare;
    compare.strings = strings;
    compare.key = key;

    void *result = bsearch(&compare,
                           entries.begin,
                           entries.size,
                           sizeof(metadata_entry_t),
                           compare_t::cmp);

    return result ? strings + ((metadata_entry_t *)result)->value : nullptr;
}

FilePointerGuard::FilePointerGuard(FilePointerPool *pool_,
                                   FILE *f_)
    : pool(pool_)
    , f(f_)
{
}

FilePointerGuard::~FilePointerGuard() {
    if(f_managed && *f_managed) {
        *f_managed = nullptr;
        f_managed = nullptr;

        pool->release(f);
    }
}

void FilePointerGuard::manage(FILE **field) {
    f_managed = field;
    *f_managed = f;
}

void FilePointerGuard::manual_release() {
    if(f_managed && *f_managed) {
        *f_managed = nullptr;
        f_managed = nullptr;

        pool->release(f);
    }
}

FilePointerPool::FilePointerPool(const std::string &path_)
    : path(path_)
{
}

FilePointerPool::~FilePointerPool() {
    for(auto f: handles) {
        fclose(f);
    }
}

FilePointerGuard FilePointerPool::acquire() {
    FILE *f = nullptr;

    while (lock.test_and_set(std::memory_order_acquire))  // acquire lock
        ; // spin

    if(!handles.empty()) {
        f = handles.back();
        handles.pop_back();
    }

    lock.clear(std::memory_order_release);               // release lock

    if(f == nullptr) {
        errif(nullptr == (f = fopen(path.c_str(), "r")),
              "Failed opening %s", path.c_str());
    }

    return FilePointerGuard(this, f);
}

void FilePointerPool::release(FILE *f) {
    while (lock.test_and_set(std::memory_order_acquire))  // acquire lock
        ; // spin

    handles.push_back(f);

    lock.clear(std::memory_order_release);               // release lock
}

PnaSequenceReader::PnaSequenceReader(FilePointerGuard fguard_,
                                     const sequence_t &sequence_,
                                     const PnaMetadata &metadata_,
                                     uint32_t flags_)
    : fguard(fguard_)
    , sequence(sequence_)
    , metadata(metadata_)
    , flags(flags_)
{
    fguard.manage(&fpna);

    seqfragments.begin = new seqfragment_t[sequence.seqfragments_count];
    errif(0 != fseeko(fpna, sequence.seqfragments_filepos, SEEK_SET),
          "Failed seeking to seqfragments");
    errif(1 != fread(seqfragments.begin, sizeof(seqfragment_t) * sequence.seqfragments_count, 1, fpna),
          "Failed reading seqfragments");
    seqfragments.next = sequence.seqfragments_count > 0 ? &seqfragments.begin[0] : nullptr;
    seqfragments.end = &seqfragments.begin[sequence.seqfragments_count];

    packedCache.buf = new unsigned char[READBUF_CAPACITY];

    this->packedByteLookup = new uint32_t[256];
    // todo: safer init
    char baseChar[] = {'A', 'C', 'G', 'T'};
    for(unsigned char b0 = 0; b0 < 4; b0++) {
        for(unsigned char b1 = 0; b1 < 4; b1++) {
            for(unsigned char b2 = 0; b2 < 4; b2++) {
                for(unsigned char b3 = 0; b3 < 4; b3++) {
                    unsigned char baseByte = (b0 << 6) | (b1 << 4) | (b2 << 2) | (b3);
                    char baseChars[] = {baseChar[b3], baseChar[b2], baseChar[b1], baseChar[b0]};
                    memcpy(packedByteLookup + baseByte, baseChars, 4);
                }
            }
        }
    }

    errif(0 != fseeko(fpna, sequence.packed_bases_filepos, SEEK_SET),
          "Failed seeking to bases");
}

PnaSequenceReader::~PnaSequenceReader() {
    delete seqfragments.begin;
    delete packedCache.buf;
    delete packedByteLookup;
}

void PnaSequenceReader::close() {
    fguard.manual_release();
}

uint64_t PnaSequenceReader::size() {
    return sequence.bases_count;
}

#define NEXT_BYTE()                                                        \
    if(packedCache.index == packedCache.len) {                                \
        errif((packedCache.bases_offset+packedCache.len) >= sequence.packed_bases_length, \
              "Attempting to read base byte when none remain!");        \
        packedCache.bases_offset += packedCache.len;                        \
        packedCache.len = min(uint64_t(READBUF_CAPACITY),                \
                              sequence.packed_bases_length - packedCache.bases_offset); \
        errif(1 != fread(packedCache.buf, packedCache.len, 1, fpna),        \
              "Failed filling read buffer.");                                \
        packedCache.index = 0;                                                \
    }                                                                        \
    packedCache.curr = packedCache.buf[packedCache.index++];

seqfragment_t *PnaSequenceReader::find_next_seqfragment(uint64_t offset) {
    struct local {
        static bool comp(const seqfragment_t &a, const seqfragment_t &b) {
            return (a.sequence_offset+a.bases_count) < b.sequence_offset;
        }
    };
    seqfragment_t searchval;
    searchval.sequence_offset = offset;
    seqfragment_t *result = lower_bound(seqfragments.begin, seqfragments.end,
                                        searchval, local::comp);

    return result == seqfragments.end ? nullptr : result;
}

void PnaSequenceReader::seek(uint64_t seekOffset) {
    if(sequence.seqfragments_count == 0) {
        seqOffset = seekOffset;
        return;
    }

    seqfragments.next = find_next_seqfragment(seekOffset);

#define RESET_CACHE(BASES_OFFSET)                        \
    packedCache.len = 0;                                \
    packedCache.index = 0;                                \
    packedCache.bases_offset = packed_bases_offset
    

    uint64_t packed_bases_offset;
    if(seqfragments.next) {
        packed_bases_offset = seqfragments.next->packed_bases_offset;

        if(seekOffset < seqfragments.next->sequence_offset) {
            // We're in an N region. Use shift state of next fragment.
            shift = seqfragments.next->shift;
        } else {
            // We're in a fragment. We must compute the new shift state
            // and the location within the packed bases.
            uint64_t relOffset = seekOffset - seqfragments.next->sequence_offset;
            shift = uint8_t( (((seqfragments.next->shift / 2) + relOffset) % 4) * 2 );
            uint64_t nfirstbyte = uint64_t((4 - (seqfragments.next->shift / 2)) % 4);
            if(relOffset >= nfirstbyte) {
                if(nfirstbyte > 0) {
                    relOffset -= nfirstbyte;
                    packed_bases_offset++;
                }
                packed_bases_offset += relOffset / 4;
            }
        }


        uint64_t packed_bases_filepos = sequence.packed_bases_filepos + packed_bases_offset;
        errif(0 != fseeko(fpna, packed_bases_filepos, SEEK_SET),
              "Failed seeking to %lu.", packed_bases_filepos);

        RESET_CACHE(packed_bases_offset);

        if(shift) {
            NEXT_BYTE();
        }
    } else {
        // No fragments left.
        packed_bases_offset = sequence.packed_bases_length;

        RESET_CACHE(packed_bases_offset);
    }

    seqOffset = seekOffset;

#undef RESET_CACHE
}

#define UNPACK(NBASES) {                                        \
        for(int i = 0; i < NBASES; i++) {                        \
            if(shift == 0) {                                        \
                NEXT_BYTE();                                        \
            }                                                        \
            switch(base_t((packedCache.curr >> shift) & 0x3)) {        \
            case A:                                                \
                *buf = 'A';                                        \
                break;                                                \
            case C:                                                \
                *buf = 'C';                                        \
                break;                                                \
            case G:                                                \
                *buf = 'G';                                        \
                break;                                                \
            case T:                                                \
                *buf = 'T';                                        \
                break;                                                \
            default:                                                \
                err("logic error");                                \
                break;                                                \
            }                                                        \
                                                                \
            shift += 2;                                                \
            buf++;                                                \
        }                                                        \
        if(shift == 8) {                                        \
            shift = 0;                                                \
        }                                                        \
    }

uint64_t PnaSequenceReader::read(char * buf, uint64_t buflen) {
    char *buf0 = buf;
    uint64_t endOffset = seqOffset + min(buflen, sequence.bases_count - seqOffset);

    while(seqOffset < endOffset) {
#define FILL_N(COUNT)                                                        \
        if(flags & IgnoreN) {                                                \
            endOffset = min(endOffset + COUNT, sequence.bases_count);        \
        } else {                                                        \
            memset(buf, 'N', COUNT);                                        \
            buf += COUNT;                                                \
        }                                                                \
        seqOffset += COUNT

        if(seqfragments.next == nullptr) {
            // The remainder of the sequence is 'N'
            size_t ncount = endOffset - seqOffset;
            FILL_N(ncount);
        } else if(seqOffset < seqfragments.next->sequence_offset) {
            // We're in an 'N' region, which is followed by a fragment.
            size_t ncount = min(endOffset - seqOffset,
                                seqfragments.next->sequence_offset - seqOffset);
            FILL_N(ncount);
        }
#undef FILL_N

        // Unpack bases while fragments remain, we're within the bounds of a fragment, and the
        // user's buffer has space.
        while(seqfragments.next && (seqOffset >= seqfragments.next->sequence_offset) && (seqOffset < endOffset)) {
            uint64_t fragment_bases_count = min(seqfragments.next->bases_count - (seqOffset - seqfragments.next->sequence_offset),
                                                endOffset - seqOffset);
            seqOffset += fragment_bases_count;

            // Unpack 1 base at a time until shift is 0
            {
                int n = (int)min(fragment_bases_count, uint64_t((4 - (shift / 2)) % 4));
                UNPACK(n);
                fragment_bases_count -= n;
            }

            // Unpack 4 bases at a time
            {
                uint32_t *intbuf = (uint32_t *)buf;
                uint64_t n = fragment_bases_count / sizeof(uint32_t);
                uint32_t *intbuf_end = intbuf + n;
                while(intbuf < intbuf_end) {
                    NEXT_BYTE();
                    *intbuf = packedByteLookup[packedCache.curr];
                    intbuf++;
                }
                buf = (char *)intbuf;
            }

            // Unpack 1 base at a time for remainder
            {
                int n = int(fragment_bases_count & 3);
                UNPACK(n);
            }

            if(seqOffset == (seqfragments.next->sequence_offset + seqfragments.next->bases_count)) {
                // We're at the end of the fragment.
                seqfragments.next++;
                if(seqfragments.next == seqfragments.end) {
                    seqfragments.next = nullptr;
                }
            }
        }
    }

    return buf - buf0;
}

const PnaMetadata PnaSequenceReader::getMetadata() {
    return metadata;
}

PnaSequenceReader::packed_read_result_t PnaSequenceReader::packed_read(uint8_t *packed_buf, uint64_t buflen) {
    errif(buflen < sequence.packed_bases_length, "Buffer too small.");

    packed_read_result_t result;

    result.bases_count = sequence.bases_count;
    result.seqfragments.begin = seqfragments.begin;
    result.seqfragments.count = sequence.seqfragments_count;
    result.packed_bases.begin = packed_buf;
    result.packed_bases.buflen = sequence.packed_bases_length;

    if(result.seqfragments.count) {
        seqfragment_t *last = seqfragments.end - 1;
        result.packed_bases.count =
            (last->packed_bases_offset * 4) + (last->shift / 2) + last->bases_count;
    }

    errif(0 != fseeko(fpna, sequence.packed_bases_filepos, SEEK_SET),
          "Failed seeking to bases");
    errif(1 != fread(packed_buf, sequence.packed_bases_length, 1, fpna),
          "Failed reading packed bases");

    return result;
}

#undef UNPACK
#undef NEXT_BYTE

PnaReader::PnaReader(const char *path_)
    : path(path_)
    , fpool(path)
{
    FILE *f;
    FilePointerGuard guard = fpool.acquire();
    guard.manage(&f);

    //
    // Read header
    //
    errif(0 != fseeko(f, 0, SEEK_SET),
          "Failed seeking to header of %s.", path.c_str());
    errif(1 != fread(&header, sizeof(header), 1, f),
          "Failed reading header of %s.", path.c_str());

    if(header.signature != PNA_FILE_SIGNATURE) {
        raise_io("PNA file signature not found.");
    }
    if(header.version != PNA_VERSION) {
        raise_io("Unsupported PNA version: %zu", size_t(header.version));
    }

    //
    // mmap strings, metadata, and sequence_t
    //
    {
        int fd = fileno(f);
        errif(fd < 0, "Failed getting fd");

        long page_size = sysconf(_SC_PAGE_SIZE);
        errif(page_size < 1, "Failed determining system page size.");

        uint64_t offset = header.string_storage.filepos;
        uint64_t end = header.sequences_filepos + (header.sequences_count * sizeof(sequence_t));
        offset = (offset / page_size) * page_size;
        mmap.length = end - offset;

        mmap.addr = ::mmap(NULL,
                           mmap.length,
                           PROT_READ,
                           MAP_SHARED,
                           fd,
                           offset);
        errif(mmap.addr == (void *)-1, "Failed mmap'ing strings");
        mmap.file_start = (uint8_t *)mmap.addr - offset;

        strings = (const char *)mmap.file_start + header.string_storage.filepos;
        sequences = (const sequence_t *)(mmap.file_start + header.sequences_filepos);
    }

    metadata = PnaMetadata(header.metadata,
                           mmap.file_start,
                           strings);
}

PnaReader::~PnaReader() {
    if(mmap.addr) {
        errif(0 != munmap(mmap.addr, mmap.length), "Failed munmap'ing");
    }
}

uint64_t PnaReader::getSequenceCount() {
    return header.sequences_count;
}

uint64_t PnaReader::getMaxSequenceFragments() {
    return header.max_seqfragments_count;
}

uint64_t PnaReader::getMaxPackedBasesLength() {
    return header.max_packed_bases_length;
}

const PnaMetadata PnaReader::getMetadata() {
    return metadata;
}

const PnaMetadata PnaReader::getSequenceMetadata(uint64_t index) {
    errif(index >= header.sequences_count,
          "Index out of bounds");

    return PnaMetadata(sequences[index].metadata,
                       mmap.file_start,
                       strings);
}

shared_ptr<PnaSequenceReader> PnaReader::openSequence(uint64_t index, uint32_t flags) {
    errif(index >= header.sequences_count,
          "Index out of bounds");

    const sequence_t &sequence = sequences[index];
    // make_shared is causing internal compiler error (gcc 4.7.3)
    return shared_ptr<PnaSequenceReader>(
        new PnaSequenceReader(fpool.acquire(),
                              sequence,
                              getSequenceMetadata(index),
                              flags));
}

uint32_t StringStorageWriter::getOffset(stringid_t id) {
    OffsetMap::iterator it = offsets.find(id);
    errif(it == offsets.end(),
          "String id not found!");

    return it->second;
}

stringid_t StringStorageWriter::getId(const char *str) {
    stringid_t result;
    String2Id::iterator it = ids.find(str);
    if(it != ids.end()) {
        result = it->second;
    } else {
        result = ids.size() + 1;
        ids[str] = result;
    }
    return result;
}

void StringStorageWriter::write(FILE *f, string_storage_t &header) {
    header.filepos = ftello(f);

    uint32_t offset = 0;
    for(auto &idpair: ids) {
        const string &str = idpair.first;
        stringid_t id = idpair.second;
        uint32_t len = str.length() + 1;

        errif((uint64_t(offset) + len) > MAX_STRING_STORAGE,
              "String storage capacity exceeded.");

        errif(1 != fwrite(str.c_str(), len, 1, f),
              "Failed writing string storage");
        offsets[id] = offset;
        offset += len;
    }

    errif((ftello(f) - header.filepos) != offset,
          "Ended at invalid location in building string storage");

    header.length = offset;
}

MetadataWriter::MetadataWriter(metadata_t &metadata_,
                               StringStorageWriter &strings_)
    : metadata(metadata_)
    , strings(strings_)
{
}

void MetadataWriter::addMetadata(const char *key, const char *value) {
    idmap[strings.getId(key)] = strings.getId(value);
}

void MetadataWriter::write(FILE *f) {
    uint32_t nentries = idmap.size();

    metadata.entries_filepos = ftello(f);
    metadata.entries_count = nentries;

    if(nentries) {
        metadata_entry_t entries[nentries];
        
        uint32_t i = 0;
        for(auto kvpair: idmap) {
            entries[i].key = strings.getOffset(kvpair.first);
            entries[i].value = strings.getOffset(kvpair.second);
            i++;
        }

        // We sort using the offset of the string into the string storage area
        // since the strings are alphabetically sorted in that region.
        struct local {
            static bool comp(const metadata_entry_t &a, const metadata_entry_t &b) {
                return a.key < b.key;
            }
        };

        sort(entries, entries + nentries, local::comp);

        errif(1 != fwrite(entries, sizeof(entries), 1, f),
              "Failed writing metadata entries");
    }
}

PnaSequenceWriter::PnaSequenceWriter(FILE *f,
                                     sequence_t &sequence_,
                                     MetadataWriter &metadata_)
    : fpna(f)
    , sequence(sequence_)
    , metadata(metadata_)
{
    // todo: don't need to normalize?
    for(int i = 0; i < 256; i++) {
        base_map[i] = N;
    }
    base_map['A'] = A;
    base_map['a'] = A;
    base_map['C'] = C;
    base_map['c'] = C;
    base_map['G'] = G;
    base_map['g'] = G;
    base_map['T'] = T;
    base_map['t'] = T;

    sequence.packed_bases_filepos = ftello(fpna);
}

PnaSequenceWriter::~PnaSequenceWriter() {
    close();
}

void PnaSequenceWriter::write(char const *buf, uint64_t buflen) {
    errif(!fpna, "Sequence closed.");

    for(uint64_t bufOffset = 0; bufOffset < buflen; bufOffset++, seqOffset++) {
        base_t base = base_map[(uint8_t)buf[bufOffset]];
        if(base != N) {
#define START_FRAGMENT()                                                \
            in_seqfragment = true;                                        \
            seqfragment.sequence_offset = seqOffset;                                \
            seqfragment.packed_bases_offset = (ftello(fpna) + packedCache.len) - sequence.packed_bases_filepos;        \
            seqfragment.shift = shift;                                        \
            seqfragment.bases_count = 1;                                        \
            
            if(!in_seqfragment) {
                START_FRAGMENT();
            } else {
                if(seqfragment.bases_count < MAX_SEQFRAGMENT_LEN) {
                    seqfragment.bases_count++;
                } else {
                    seqfragments.push_back(seqfragment);
                    START_FRAGMENT();
                }
            }
#undef START_FRAGMENT

            packedByte |= base << shift;

            shift += 2;
            if(shift == 8) {
                addPackedByte();
                shift = 0;
                packedByte = 0;
            }

        } else {
            if(in_seqfragment) {
                in_seqfragment = false;
                seqfragments.push_back(seqfragment);
            }
        }
    }
}

void PnaSequenceWriter::addMetadata(const char *key, const char *value) {
    errif(!fpna, "Sequence closed.");

    metadata.addMetadata(key, value);
}

uint64_t PnaSequenceWriter::getBasePairCount() {
    errif(!fpna, "Sequence closed");
    return seqOffset;
}

uint64_t PnaSequenceWriter::getByteCount() {
    errif(!fpna, "Sequence closed");
    return
        sizeof(sequence_t)
        + sizeof(seqfragment_t) * seqfragments.size()
        + ftello(fpna) - sequence.packed_bases_filepos;
}

void PnaSequenceWriter::close() {
    if(!fpna)
        return;

    if(shift != 0) {
        addPackedByte();
    }
    flushCache();

    sequence.bases_count = seqOffset;
    sequence.packed_bases_length = ftello(fpna) - sequence.packed_bases_filepos;

    if(in_seqfragment) {
        seqfragment.bases_count = seqOffset - seqfragment.sequence_offset;
        seqfragments.push_back(seqfragment);
    }

    sequence.seqfragments_filepos = ftello(fpna);
    sequence.seqfragments_count = seqfragments.size();

    for(seqfragment_t seqfragment: seqfragments) {
        errif(1 != fwrite(&seqfragment, sizeof(seqfragment_t), 1, fpna),
              "Failed writing seqfragment");
    }

    fpna = nullptr;
}

void PnaSequenceWriter::addPackedByte() {
    if(packedCache.len == WRITEBUF_CAPACITY) {
        flushCache();
    }
    packedCache.buf[packedCache.len++] = packedByte;
}

void PnaSequenceWriter::flushCache() {
    if(packedCache.len) {
        errif(1 != fwrite(packedCache.buf, packedCache.len, 1, fpna),
              "Failed writing packed cache!!!\n");
        packedCache.len = 0;
    }
}

PnaWriter::PnaWriter(const char *path) // todo: const string &
    : metadata(header.metadata, strings)
{
    f = fopen(path, "w");
    memset(&header, 0, sizeof(header));

    header.signature = PNA_FILE_SIGNATURE;
    header.version = PNA_VERSION;

    errif(0 != fseeko(f, sizeof(header), SEEK_SET),
          "Failed seeking past header");
    header.sequences_filepos = ftello(f);
}

PnaWriter::~PnaWriter() {
    close();
}

shared_ptr<PnaSequenceWriter> PnaWriter::createSequence() {
    closeSequence();

    sequence_t *sequence = new sequence_t;
    MetadataWriter *metadata = new MetadataWriter(sequence->metadata, strings);

    sequences.push_back(sequence);
    sequences_metadata.push_back(metadata);

    shared_ptr<PnaSequenceWriter> writer(new PnaSequenceWriter(f, *sequence, *metadata));
    activeSequenceWriter = writer;
    return writer;
}

void PnaWriter::closeSequence() {
    if(auto prev = activeSequenceWriter.lock()) {
        prev->close();
        activeSequenceWriter.reset();
    }
}

void PnaWriter::addMetadata(const char *key, const char *value) {
    errif(!f, "File closed!");

    metadata.addMetadata(key, value);
}

void PnaWriter::close() {
    if(!f)
        return;

    closeSequence();

    //
    // Write strings storage
    //
    strings.write(f, header.string_storage);

    //
    // Write metadata entries
    //
    for(auto metadata: sequences_metadata) {
        metadata->write(f);
        delete metadata;
    }
    sequences_metadata.clear();

    metadata.write(f);

    //
    // Write sequence_t set
    //
    header.sequences_filepos = ftello(f);
    header.sequences_count = sequences.size();

    for(auto sequence: sequences) {
        errif(1 != fwrite(sequence, sizeof(sequence_t), 1, f),
              "Failed writing sequence");
        header.max_seqfragments_count = max(header.max_seqfragments_count,
                                            sequence->seqfragments_count);
        header.max_packed_bases_length = max(header.max_packed_bases_length,
                                             sequence->packed_bases_length);
        delete sequence;
    }
    sequences.clear();

    //
    // Write header
    //
    errif(0 != fseeko(f, 0, SEEK_SET),
          "Failed seeking to header");
    errif(1 != fwrite(&header, sizeof(header), 1, f),
          "Failed writing header");

    //
    // Done
    //
    fclose(f);
    f = nullptr;
}
