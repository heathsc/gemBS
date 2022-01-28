/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_file.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_FILE_H_
#define GT_INPUT_FILE_H_

#include "gt_essentials.h"
//#include "gt_attributes.h"
//#include "gt_sam_attributes.h"

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif
#ifdef HAVE_BZLIB
#include <bzlib.h>
#endif

// Codes gt_status
#define GT_INPUT_FILE_OK 0
#define GT_INPUT_FILE_CLOSE_ERR 1
#define GT_INPUT_FILE_EOF 0
#define GT_INPUT_FILE_LINE_READ 1

/*
 * Checkers
 */
#define GT_INPUT_FILE_CHECK(input_file) \
  GT_NULL_CHECK(input_file); \
  GT_NULL_CHECK(input_file->file_name)

/*
 * GT Input file
 */
typedef enum { FASTA, FASTQ, MAP, SAM, GENERIC, FILE_FORMAT_UNKNOWN } gt_file_format;
typedef enum { STREAM, REGULAR_FILE, MAPPED_FILE, GZIPPED_FILE, BZIPPED_FILE, XZIPPED_FILE } gt_file_type;
typedef struct {
  /* Input file */
  char* file_name;
  gt_file_type file_type;
  FILE* file;
  int fildes;
  bool eof;
  uint64_t file_size;
  pthread_mutex_t input_mutex;
  /* Auxiliary Buffer (for synch purposes) */
  uint8_t* file_buffer;
  uint64_t buffer_size;
  uint64_t buffer_begin;
  uint64_t buffer_pos;
  uint64_t global_pos;
  uint64_t processed_lines;
  /* ID generator */
  uint64_t processed_id;
} gt_input_file;

/*
 * Basic I/O functions
 */
gt_input_file* gt_input_stream_general_open(FILE* stream,gt_file_format format);
gt_input_file* gt_input_file_general_open(char* const file_name,const bool mmap_file,gt_file_format format);
gt_status gt_input_file_close(gt_input_file* const input_file);

#define gt_input_stream_open(stream) gt_input_stream_general_open(stream,FILE_FORMAT_UNKNOWN)
#define gt_input_stream_map_open(stream) gt_input_stream_general_open(stream,MAP)
#define gt_input_stream_fasta_open(stream) gt_input_stream_general_open(stream,FASTA)
#define gt_input_stream_sam_open(stream) gt_input_stream_general_open(stream,SAM)

#define gt_input_file_open(file_name,mmap_file) gt_input_file_general_open(file_name,mmap_file,FILE_FORMAT_UNKNOWN)
#define gt_input_file_map_open(file_name,mmap_file) gt_input_file_general_open(file_name,mmap_file,MAP)
#define gt_input_file_fasta_open(file_name,mmap_file) gt_input_file_general_open(file_name,mmap_file,FASTA)
#define gt_input_file_sam_open(file_name,mmap_file) gt_input_file_general_open(file_name,mmap_file,SAM)
#define gt_input_file_generic_open(file_name,mmap_file) gt_input_file_general_open(file_name,mmap_file,GENERIC)

/*
 * Accessors (Mutex,ID,...) functions
 */
GT_INLINE void gt_input_file_lock(gt_input_file* const input_file);
GT_INLINE void gt_input_file_unlock(gt_input_file* const input_file);
GT_INLINE uint64_t gt_input_file_next_id(gt_input_file* const input_file);

/*
 * Basic line functions
 */
GT_INLINE size_t gt_input_file_dump_to_buffer(gt_input_file* const input_file,gt_vector* const buffer_dst);
GT_INLINE size_t gt_input_file_fill_buffer(gt_input_file* const input_file);
GT_INLINE size_t gt_input_file_next_line(gt_input_file* const input_file,gt_vector* const buffer_dst);

GT_INLINE size_t gt_input_file_next_record(
    gt_input_file* const input_file,gt_vector* const buffer_dst,gt_string* const first_field,
    uint64_t* const num_spaces,uint64_t* const num_tabs);
GT_INLINE bool gt_input_file_next_record_cmp_first_field(gt_input_file* const input_file,gt_string* const first_field);

/*
 * Line Readers (thread-unsafe, must call mutex functions before)
 */
GT_INLINE uint64_t gt_input_file_add_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines);
GT_INLINE uint64_t gt_input_file_get_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines);

/*
 * Processing Macros (direct parsing from input file)
 */
#define GT_INPUT_FILE_CHECK_BUFFER(input_file) \
  if (gt_expect_false(input_file->buffer_pos >= input_file->buffer_size)) { \
    gt_input_file_fill_buffer(input_file); \
  }
#define GT_INPUT_FILE_CHECK_BUFFER__DUMP(input_file,buffer_dst) \
  if (gt_expect_false(input_file->buffer_pos >= input_file->buffer_size)) { \
    if (gt_expect_true(buffer_dst!=NULL)) gt_input_file_dump_to_buffer(input_file,buffer_dst); \
    gt_input_file_fill_buffer(input_file); \
  }

#define GT_INPUT_FILE_NEXT_CHAR(input_file) \
  ++input_file->buffer_pos; \
  GT_INPUT_FILE_CHECK_BUFFER(input_file)
#define GT_INPUT_FILE_NEXT_CHAR__DUMP(input_file,buffer_dst) \
  ++input_file->buffer_pos; \
  GT_INPUT_FILE_CHECK_BUFFER__DUMP(input_file,buffer_dst)

#define GT_INPUT_FILE_SKIP_EOL(input_file) \
  if (!input_file->eof) { \
    GT_INPUT_FILE_NEXT_CHAR(input_file); /* Skip EOF/DOS_EOL */  \
    if (gt_expect_false(!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)==EOL)) { \
      GT_INPUT_FILE_NEXT_CHAR(input_file); /* Skip DOS_EOL */  \
    } \
  }
#define GT_INPUT_FILE_HANDLE_EOL(input_file,buffer_dst) \
  if (!input_file->eof) { \
    GT_INPUT_FILE_NEXT_CHAR__DUMP(input_file,buffer_dst); /* Skip EOF/DOS_EOL */  \
    if (gt_expect_false(!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)==EOL)) { \
      ++input_file->buffer_pos; \
      if (gt_expect_true(buffer_dst!=NULL)) { \
        gt_input_file_dump_to_buffer(input_file,buffer_dst); \
        gt_vector_dec_used(buffer_dst); \
        *gt_vector_get_last_elm(buffer_dst,char)=EOL; \
      } \
      if (gt_expect_false(input_file->buffer_pos >= input_file->buffer_size)) { \
        gt_input_file_fill_buffer(input_file); \
      } \
    } \
  }

#define GT_INPUT_FILE_CURRENT_CHAR(input_file) input_file->file_buffer[input_file->buffer_pos]


#endif /* GT_INPUT_FILE_H_ */
