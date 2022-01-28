/*
 * PROJECT: GEM-Tools library
 * FILE: gt_fm.h
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   // TODO
 */

#ifndef GT_FILE_MANAGEMENT_H_
#define GT_FILE_MANAGEMENT_H_

#include "gt_commons.h"
#include "gt_error.h"
#include "gt_vector.h"
#include "gt_string.h"

#include <unistd.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/*
 * I/O Constants/Values
 */
extern int gt_fm_oflags[3];

typedef enum { GT_FILE_READ_ONLY=0, GT_FILE_WRITE_ONLY=1, GT_FILE_READ_WRITE=2 } gt_fm_mode;
typedef struct {
  /* File */
  int fd;                 /* File descriptor */
  FILE* file;             /* FILE */
  /* Attributes */
  char *file_name;        /* File name */
  gt_fm_mode mode;        /* File mode {R,W,R/W} */
  /* Locator */
  uint64_t byte_position; /* Current byte position */
} gt_fm;

/*
 * Setup
 */
GT_INLINE gt_fm* gt_fm_open_file(char* const file_name,const gt_fm_mode mode);
GT_INLINE gt_fm* gt_fm_open_fd(const int fd,const gt_fm_mode mode);
GT_INLINE gt_fm* gt_fm_open_FILE(FILE* const stream,const gt_fm_mode mode);

GT_INLINE void gt_fm_close(gt_fm* const fm);

/*
 * Accessors
 */
GT_INLINE uint64_t gt_fm_get_current_position(gt_fm* const fm);
GT_INLINE bool gt_fm_eof(gt_fm* const fm);
GT_INLINE gt_fm_mode gt_fm_get_mode(gt_fm* const fm);
GT_INLINE void gt_fm_set_mode(gt_fm* const fm,const gt_fm_mode mode);

/*
 * Bulk read of file
 */
GT_INLINE void gt_fm_bulk_read_fd(const int fd,void* const dst,const uint64_t size);
GT_INLINE void gt_fm_bulk_read_file_parallel(
    char* const file_name,void* const dst,const uint64_t offset,const uint64_t size,const uint64_t num_threads);
GT_INLINE void gt_fm_bulk_read_file(char* const file_name,void* const dst,const uint64_t offset,const uint64_t size);

/*
 * Seek
 */
GT_INLINE void gt_fm_seek(gt_fm* const fm,const uint64_t byte_position);
GT_INLINE void gt_fm_skip_forward(gt_fm* const fm,const uint64_t num_bytes);
GT_INLINE void gt_fm_skip_backward(gt_fm* const fm,const uint64_t num_bytes);
GT_INLINE void gt_fm_skip_uint64(gt_fm* const fm);
GT_INLINE void gt_fm_skip_uint32(gt_fm* const fm);
GT_INLINE void gt_fm_skip_uint16(gt_fm* const fm);
GT_INLINE void gt_fm_skip_uint8(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align(gt_fm* const fm,const uint64_t num_bytes);
GT_INLINE void gt_fm_skip_align_16(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_32(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_64(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_128(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_512(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_1024(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_4KB(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_mempage(gt_fm* const fm);

/*
 * Read
 */
#define gt_fm_read(fm,var) gt_fm_copy_mem(fm,&var,sizeof(var))
GT_INLINE uint64_t gt_fm_read_uint64(gt_fm* const fm);
GT_INLINE uint32_t gt_fm_read_uint32(gt_fm* const fm);
GT_INLINE uint16_t gt_fm_read_uint16(gt_fm* const fm);
GT_INLINE uint8_t gt_fm_read_uint8(gt_fm* const fm);
GT_INLINE void* gt_fm_read_mem(gt_fm* const fm,const uint64_t num_bytes);
GT_INLINE void gt_fm_copy_mem(gt_fm* const fm,void* const dst,const uint64_t num_bytes);
GT_INLINE void gt_fm_copy_mem_parallel(gt_fm* const fm,void* const dst,const uint64_t num_bytes,const uint64_t num_threads);

/*
 * Write
 */
#define gt_fm_write(fm,var) gt_fm_write_mem(fm,&var,sizeof(var))
GT_INLINE void gt_fm_write_uint64(gt_fm* const fm,const uint64_t data);
GT_INLINE void gt_fm_write_uint32(gt_fm* const fm,const uint32_t data);
GT_INLINE void gt_fm_write_uint16(gt_fm* const fm,const uint16_t data);
GT_INLINE void gt_fm_write_uint8(gt_fm* const fm,const uint8_t data);
GT_INLINE void gt_fm_write_mem(gt_fm* const fm,void* const src,const uint64_t num_bytes);

#endif /* GT_FILE_MANAGEMENT_H_ */

