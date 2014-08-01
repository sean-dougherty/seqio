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
  Options passed to seqio_get_sequences().
 */
typedef struct {
    seqio_base_transform base_transform;
} seqio_sequence_options;

/*!
  The possible return values for seqio_get_status()
*/
typedef enum {
    SEQIO_SUCCESS = 0, 
    SEQIO_ERR_INVALID_PARAMETER = 1, 
    SEQIO_ERR_FILE_NOT_FOUND = 2,
    SEQIO_ERR_IO = 3,
    SEQIO_ERR_KEY_NOT_FOUND = 4,
    SEQIO_ERR_OUT_OF_MEMORY = 5
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
  to seqio_dispose_sequences().
*/
typedef void *seqio_sequence_iterator;

/*!
  Handle to an object representing a sequence, which can be used to obtain metadata
  (e.g. sequence name, comment) and sequence data.
*/
typedef void *seqio_sequence;

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
  Provides reasonable default options for seqio_get_sequences():
  - base_transform: SEQIO_BASE_TRANSFORM_NONE
*/
extern seqio_sequence_options const SEQIO_DEFAULT_SEQUENCE_OPTIONS;

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
  seqio_get_sequences("./test.fa.gz",
                      SEQIO_DEFAULT_SEQUENCE_OPTIONS,
                      &iterator);

  // Get the first sequence.
  seqio_next_sequence(iterator, &sequence);

  // Print out its name.
  seqio_get_value(sequence, SEQIO_KEY_NAME, &name);
  
  \endcode
*/
seqio_status seqio_get_sequences(char const *path,
                                 seqio_sequence_options options,
                                 seqio_sequence_iterator *iterator);

/*!
  Dispose the iterator and any sequences that were obtained from it.
   
  \param [in,out] iterator The iterator to be disposed. Value is set to SEQIO_NIL_SEQUENCE_ITERATOR.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.

  \note SEQIO_NIL_SEQUENCE_ITERATOR is accepted and is a no-op.
*/
seqio_status seqio_dispose_sequences(seqio_sequence_iterator *iterator);

/*!
  Specifies if seqio_next_sequence() will return another sequence.

  \param [in] iterator The iterator being tested.
  \param [out] has_next On return, specifies whether another sequence exists.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.
*/
    seqio_status seqio_has_next_sequence(seqio_sequence_iterator iterator,
                                         seqio_bool *has_next);

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
  seqio_bool has_next;

  char *buffer = NULL;
  uint32_t buffer_length;

  // ... create iterator ...

  while( (seqio_has_next_sequence(iterator, &has_next) == SEQIO_SUCCESS)
         && has_next) {
    uint32_t sequence_length;

    seqio_next_sequence(iterator, &sequence);

    // By reading the sequence, the file pointer will be at the end of the sequence,
    // allowing the iterator to efficiently find the next sequence in a format like FASTA.
    seqio_read_all(sequence, &sequence_length, &buffer, &buffer_length);

    // ... do stuff with sequence data ...

    seqio_dispose(sequence);
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
    seqio_status seqio_dispose(seqio_sequence *sequence);

/*!
  Get the number of metadata keys for the sequence.

  \param [in] sequence The sequence whose metadata is being queried.
  \param [out] count Upon return, the number of keys in the sequence's metadata.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.
 */
    seqio_status seqio_get_key_count(seqio_sequence sequence,
                                     uint32_t *count);

/*!
  Get the name of the key at specified index.

  \param [in] sequence The sequence whose metadata is being queried.
  \param [in] key_index Index of key being queried. Must be less than key count.
  \param [out] key Upon return, the name of the key. Client code should not free
               this.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.
*/
    seqio_status seqio_get_key(seqio_sequence sequence,
                               uint32_t key_index,
                               char const **key);
  
/*!
  Get a metadata value for the sequence.

  \param [in] sequence The sequence whose metadata is being queried.
  \param [in] key Name of key being queried.
  \param [out] value Upon return, value for the key. Client code should not free this.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.
 */
    seqio_status seqio_get_value(seqio_sequence sequence,
                                 char const *key,
                                 char const **value);

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
  Read entirety of sequence, allocating buffer on client's behalf.

  \param [in] sequence The sequence to be read.
  \param [out] sequence_length Upon return, number of bases read.
  \param [in,out] buffer Destination for bases. Client should not allocate the memory
                         pointed to by *buffer. *buffer should initially be NULL.
  \param [in,out] buffer_length Upon return, size of buffer in bytes. If *buffer parameter
                                is NULL, this value is ignored.

  \return SEQIO_SUCCESS if successful, otherwise SEQIO_ERR_*.

  \b Example
  \code
  seqio_sequence_iterator iterator;
  seqio_sequence sequence;
  char *buffer = NULL;
  uint32_t buffer_length;
  uint32_t sequence_length;

  // ... obtain sequence iterator ...

  sequence = seqio_next_sequence(iterator, &sequence);
  seqio_read_all(sequence,
                 &sequence_length,
                 &buffer,
                 &buffer_length);

  //
  
  \endcode
 */
    seqio_status seqio_read_all(seqio_sequence sequence,
                                uint32_t *sequence_length,
                                char **buffer,
                                uint32_t *buffer_length);

/*!
  Dispose a buffer allocated by the internal implementation (e.g. seqio_read_all()).

  \param [in,out] buffer The buffer to be disposed. *buffer can be NULL. Upon return,
                         will be NULL.

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
