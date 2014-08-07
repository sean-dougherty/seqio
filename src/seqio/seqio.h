#pragma once

#include <stdint.h>

/*!
  A standard key for obtaining the name of a sequence via seqio_get(). This key will be
  valid for any file format.
*/
#define SEQIO_KEY_NAME "seqio.name"

/*!
  A standard key for obtaining a FASTA-style comment for a sequence via seqio_get().
  A sequence in a PNA file will typically have a FASTA comment associated with it.
 */
#define SEQIO_KEY_COMMENT "seqio.comment"

/*!
  Value of a null seqio_sequence.
 */
#define SEQIO_NIL_SEQUENCE ((seqio_sequence *)0)

/*!
  Value of a null seqio_sequence_iterator.
 */
#define SEQIO_NIL_SEQUENCE_ITERATOR ((seqio_sequence_iterator *)0)

typedef enum {
    SEQIO_FILE_FORMAT_DEDUCE,
    SEQIO_FILE_FORMAT_FASTA,
    SEQIO_FILE_FORMAT_FASTA_GZIP,
    SEQIO_FILE_FORMAT_PNA
} seqio_file_format;

/*!
  Specifies transformation that should be applied to bases prior to being placed in
  client read buffer.
*/
typedef enum {
    /*! Don't transform bases. */
    SEQIO_BASE_TRANSFORM_NONE,
    /*! Make all bases uppercase (GACT) and make anything other than GACT an N. */
    SEQIO_BASE_TRANSFORM_CAPS_GATCN
} seqio_base_transform;

/*!
  Options passed to seqio_create_sequence_iterator().
 */
typedef struct {
    seqio_file_format file_format;
    seqio_base_transform base_transform;
} seqio_sequence_options;

typedef struct {
    seqio_file_format file_format;
} seqio_writer_options;

/*!
  The possible return values for seqio_get_status()
*/
typedef enum {
    SEQIO_SUCCESS = 0, 
    SEQIO_ERR_INVALID_PARAMETER = 1, 
    SEQIO_ERR_INVALID_STATE = 2, 
    SEQIO_ERR_FILE_NOT_FOUND = 3,
    SEQIO_ERR_IO = 4,
    SEQIO_ERR_KEY_NOT_FOUND = 5,
    SEQIO_ERR_OUT_OF_MEMORY = 6
} seqio_status;

/*!
  Boolean type that is either SEQIO_TRUE or SEQIO_FALSE
 */
typedef uint8_t seqio_bool;
/*!
  True value for type seqio_bool.
 */
#define SEQIO_TRUE 1
/*!
  False value for type seqio_bool.
 */
#define SEQIO_FALSE 0

/*!
  Handle to an object capable of iterating over a set of sequences. For example, an iterator
  for a FASTA file containing 5 sequences will return 5 seqio_sequence handles.

  Resources associated with the iterator (e.g. file handle) will not be released until a call
  to seqio_dispose_sequence_iterator().
*/
typedef struct __seqio_sequence_iterator *seqio_sequence_iterator;

typedef struct __seqio_dictionary *seqio_dictionary;
typedef struct __seqio_dictionary const *seqio_const_dictionary;

/*!
  Handle to an object representing a sequence, which can be used to obtain metadata
  (e.g. sequence name, comment) and sequence data.
*/
typedef struct __seqio_sequence *seqio_sequence;

typedef struct __seqio_writer *seqio_writer;

/*!
  Information passed to error handler.
*/
typedef struct {
    seqio_status status;
    char const *message;
} seqio_err_info;

/*!
  Signature of error handler.
*/
typedef seqio_status (*seqio_err_handler)(seqio_err_info err_info);

#ifdef __cplusplus
extern "C" {
#endif

/*!
  Provides reasonable default options for seqio_create_sequence_iterator():
  - base_transform: SEQIO_BASE_TRANSFORM_NONE
*/
extern seqio_sequence_options const SEQIO_DEFAULT_SEQUENCE_OPTIONS;
extern seqio_writer_options const SEQIO_DEFAULT_WRITER_OPTIONS;

extern seqio_err_handler const SEQIO_ERR_HANDLER_ABORT;
extern seqio_err_handler const SEQIO_ERR_HANDLER_RETURN;

/*!
  Open a file for reading one or more sequences. The file must be one of the supported formats:
  - FASTA (may be gzipped)
  - PNA

  The set of sequences will be lazily populated for inherently sequential formats like FASTA.

  \param [in] path Location of file in filesystem.
  \param [in] options Options used in processing the sequences.
  \param [out] iterator Iterator for all sequences found in file.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.

  \attention File formats like FASTA are inherently sequential, and the seqio
  implementation must therefore read a sequence's entire contents in order to
  find the next sequence. For this reason, it is most efficient to read the
  contents of sequences between calls to sequence_next().

  \b Example
  \code
  seqio_sequence_iterator iterator;
  seqio_sequence sequence;
  const char *name;

  // Open the file.
  seqio_create_sequence_iterator("./test.fa.gz",
                        SEQIO_DEFAULT_SEQUENCE_OPTIONS,
                        &iterator);

  // Get the first sequence.
  seqio_next_sequence(iterator, &sequence);

  // Print out its name.
  seqio_get_value(sequence, SEQIO_KEY_NAME, &name);
  
  \endcode
*/
seqio_status seqio_create_sequence_iterator(char const *path,
                                            seqio_sequence_options options,
                                            seqio_sequence_iterator *iterator);

/*!
  Dispose the iterator. Any sequences that were obtained from it that have not
  been disposed are still valid.
   
  \param [in,out] iterator The iterator to be disposed. Value is set to SEQIO_NIL_SEQUENCE_ITERATOR.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.

  \note SEQIO_NIL_SEQUENCE_ITERATOR is accepted and is a no-op.
*/
seqio_status seqio_dispose_sequence_iterator(seqio_sequence_iterator *iterator);

/*!
  Get the next sequence from the iterator. The result is placed in the sequence
  parameter.

  \param [in] iterator Sequence iterator.
  \param [out] sequence Contains next sequence upon return.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.

  \attention File formats like FASTA are inherently sequential, and the seqio
  implementation must therefore read a sequence's entire contents in order to
  find the next sequence. For this reason, it is most efficient to read the
  contents of sequences between calls to sequence_next(). See example code.

  \b Example
  \code
  seqio_sequence_iterator iterator;
  seqio_sequence sequence;

  char *buffer = NULL;
  uint32_t buffer_length;

  // ... create iterator ...

  while( (seqio_next_sequence(iterator, &sequence) == SEQIO_SUCCESS)
         && iterator) {
    uint32_t sequence_length;

    // By reading the sequence, the file pointer will be at the end of the sequence,
    // allowing the iterator to efficiently find the next sequence in a format like FASTA.
    seqio_read_all(sequence, &sequence_length, &buffer, &buffer_length);

    // ... do stuff with sequence data ...

    seqio_dispose_sequence(sequence);
  }

  // ... dispose iterator and buffer ...

  \endcode
 */
    seqio_status seqio_next_sequence(seqio_sequence_iterator iterator,
                                     seqio_sequence *sequence);

/*!
  Disposes resources associated with a sequence (e.g. read cache).

  \param [in,out] sequence The sequence to be disposed. Value is set to SEQIO_NIL_SEQUENCE.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.

  \note SEQIO_NIL_SEQUENCE is accepted and is a no-op.
 */
    seqio_status seqio_dispose_sequence(seqio_sequence *sequence);

    seqio_status seqio_get_metadata(seqio_sequence sequence,
                                    seqio_const_dictionary *dict);

/*!
  Read bases from sequence.

  \param [in] sequence The sequence to be read.
  \param [in] buffer Destination for bases.
  \param [in] buffer_length Size of buffer in bytes.
  \param [out] read_length Number of bases read.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.
 */
    seqio_status seqio_read(seqio_sequence sequence,
                            char *buffer,
                            uint32_t buffer_length,
                            uint32_t *read_length);

/*!
  Read entirety of sequence, allocating buffer on client's behalf. If some portion of the
  sequence has already been read via seqio_read(), that portion will not be included in the
  result of read_all(). The result is guaranteed to be null-terminated; the null is not
  included in read_length.

  \param [in] sequence The sequence to be read.
  \param [in,out] buffer Destination for bases. Client should not allocate the memory
                         pointed to by *buffer. *buffer should initially be NULL.
  \param [in,out] buffer_length Upon return, size of buffer in bytes. If *buffer parameter
                                is NULL, this value is ignored.
  \param [out] read_length Upon return, number of bases read.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.

  \b Example
  \code
  seqio_sequence_iterator iterator;
  seqio_sequence sequence;
  char *buffer = NULL;
  uint32_t buffer_length;
  uint32_t read_length;

  // ... obtain sequence iterator ...

  sequence = seqio_next_sequence(iterator, &sequence);
  seqio_read_all(sequence,
                 &read_length,
                 &buffer,
                 &buffer_length);

  //
  
  \endcode
 */
    seqio_status seqio_read_all(seqio_sequence sequence,
                                char **buffer,
                                uint32_t *buffer_length,
                                uint32_t *read_length);

    seqio_status seqio_create_writer(char const *path,
                                     seqio_writer_options options,
                                     seqio_writer *writer);

    seqio_status seqio_dispose_writer(seqio_writer *writer);

    seqio_status seqio_create_sequence(seqio_writer writer);

    seqio_status seqio_add_metadata(seqio_writer writer,
                                    char const *key,
                                    char const *value);

    seqio_status seqio_write(seqio_writer writer,
                             char const *buffer,
                             uint32_t length);

/*
    seqio_status seqio_has_key(seqio_const_dictionary dict,
                               seqio_bool *has_key);
*/

    seqio_status seqio_get_key_count(seqio_const_dictionary dict,
                                     uint32_t *count);

    seqio_status seqio_get_key(seqio_const_dictionary dict,
                               uint32_t key_index,
                               char const **key);
  
    seqio_status seqio_get_value(seqio_const_dictionary dict,
                                 char const *key,
                                 char const **value);

    seqio_status seqio_set_value(seqio_dictionary dict,
                                 char const *key,
                                 char const *value);

/*!
  Dispose a buffer allocated by the internal implementation (e.g. seqio_read_all()).

  \param [in,out] buffer The buffer to be disposed. *buffer can be NULL. Upon return,
                         *buffer will be NULL.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.
*/  
    seqio_status seqio_dispose_buffer(char **buffer);

/*!
  Set the error handler.
*/
    seqio_status seqio_set_err_handler(seqio_err_handler err_handler);

#ifdef __cplusplus
} // extern "C"
#endif
