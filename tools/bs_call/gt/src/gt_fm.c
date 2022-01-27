/*
 * PROJECT: GEM-Tools library
 * FILE: gt_fm.c
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
 *   1.- madvise() / readahead() / posix_fadvise()
 */

#include "gt_fm.h"

#ifndef O_NOATIME
  #define O_NOATIME 0
#endif
/*
 * I/O Constants/Values
 */
#define GT_FM_BULK_COPY_BLOCK_SIZE (1<<30) /* 1GB */
// FIXME: O_NOATIME not allowed is only read rights are guaranteed
int gt_fm_oflags[3] = { O_RDONLY, O_WRONLY|O_CREAT|O_NOATIME, O_RDWR|O_CREAT|O_NOATIME };

/*
 * Setup
 */
GT_INLINE gt_fm* gt_fm_open_file(char* const file_name,const gt_fm_mode mode) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE gt_fm* gt_fm_open_fd(const int fd,const gt_fm_mode mode) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE gt_fm* gt_fm_open_FILE(FILE* const stream,const gt_fm_mode mode) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}

GT_INLINE void gt_fm_close(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}

/*
 * Accesors
 */
GT_INLINE uint64_t gt_fm_get_current_position(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE bool gt_fm_eof(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE gt_fm_mode gt_fm_get_mode(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_set_mode(gt_fm* const fm,const gt_fm_mode mode) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}

/*
 * Bulk read of file
 */
GT_INLINE void gt_fm_bulk_read_fd(const int fd,void* const dst,const uint64_t size) {
  GT_NULL_CHECK(dst);
  uint64_t bytes_written = 0;
  while (bytes_written < size) {
    const uint64_t bytes_pending = size-bytes_written;
    const uint64_t chunk_size = (bytes_pending<GT_FM_BULK_COPY_BLOCK_SIZE) ? bytes_pending : GT_FM_BULK_COPY_BLOCK_SIZE;
    // Copy chunk
    gt_cond_fatal_error__perror(read(fd,dst+bytes_written,chunk_size)!=chunk_size,FILE_READ_FD);
    bytes_written += chunk_size;
  }
}
GT_INLINE void gt_fm_bulk_read_file_parallel(
    char* const file_name,void* const dst,const uint64_t offset,const uint64_t size,const uint64_t num_threads) {
  GT_NULL_CHECK(file_name);
  GT_NULL_CHECK(dst);
  // Retrieve input file info
  struct stat stat_info;
  gt_cond_fatal_error__perror(stat(file_name,&stat_info)==-1,FILE_STAT,file_name);
  // Calculate size of each chunk
  const uint64_t size_to_read = (size==0) ? stat_info.st_size : size;
  const uint64_t chunk_size = size_to_read/num_threads; // stat_info.st_size > num_threads
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(num_threads)
#endif
  {
    // Calculate offsets
#ifdef HAVE_OPENMP
    const uint64_t tid = omp_get_thread_num();
#else
    const uint64_t tid = 0;
#endif
    const uint64_t thread_mem_offset = tid*chunk_size;
    const uint64_t thread_file_offset = offset + thread_mem_offset;
    const uint64_t thread_size = (tid < (num_threads-1)) ? chunk_size : stat_info.st_size-chunk_size*tid;
    // Open file descriptor
    int fd = open(file_name,O_RDONLY,S_IRUSR);
    if(fd == -1){
    	gt_fatal_error(FILE_OPEN,file_name);
    }
    if(lseek(fd,thread_file_offset,SEEK_SET)==-1){
    	gt_fatal_error(FILE_SEEK,file_name,thread_file_offset); // Seek
    }
    // Copy file chunk
    gt_fm_bulk_read_fd(fd,dst+thread_mem_offset,thread_size);
  }
}
GT_INLINE void gt_fm_bulk_read_file(char* const file_name,void* const dst,const uint64_t offset,const uint64_t size) {
  GT_NULL_CHECK(file_name);
  GT_NULL_CHECK(dst);
  // Retrieve input file info
  struct stat stat_info;
  gt_cond_fatal_error__perror(stat(file_name,&stat_info)==-1,FILE_STAT,file_name);
  // Open file descriptor
  int fd = open(file_name,O_RDONLY,S_IRUSR);
  gt_cond_fatal_error__perror(fd==-1,FILE_OPEN,file_name);
  if (offset > 0) {
    gt_cond_fatal_error__perror(lseek(fd,offset,SEEK_SET)==-1,FILE_SEEK,file_name,offset); // Seek
  }
  // Read file
  gt_fm_bulk_read_fd(fd,dst,(size==0) ? stat_info.st_size-offset : size);
}

/*
 * Seek
 */
GT_INLINE void gt_fm_seek(gt_fm* const fm,const uint64_t byte_position) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_forward(gt_fm* const fm,const uint64_t num_bytes) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_backward(gt_fm* const fm,const uint64_t num_bytes);
GT_INLINE void gt_fm_skip_uint64(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_uint32(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_uint16(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_uint8(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_align(gt_fm* const fm,const uint64_t num_bytes) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_align_16(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_align_32(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_align_64(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_align_128(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_align_512(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_align_1024(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_align_4KB(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_skip_align_mempage(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}

/*
 * Read
 */
GT_INLINE uint64_t gt_fm_read_uint64(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE uint32_t gt_fm_read_uint32(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE uint16_t gt_fm_read_uint16(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE uint8_t gt_fm_read_uint8(gt_fm* const fm) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void* gt_fm_read_mem(gt_fm* const fm,const uint64_t num_bytes) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_copy_mem(gt_fm* const fm,void* const dst,const uint64_t num_bytes) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_copy_mem_parallel(gt_fm* const fm,void* const dst,const uint64_t num_bytes,const uint64_t num_threads) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}

/*
 * Write
 */
GT_INLINE void gt_fm_write_uint64(gt_fm* const fm,const uint64_t data) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_write_uint32(gt_fm* const fm,const uint32_t data) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_write_uint16(gt_fm* const fm,const uint16_t data) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_write_uint8(gt_fm* const fm,const uint8_t data) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
GT_INLINE void gt_fm_write_mem(gt_fm* const fm,void* const src,const uint64_t num_bytes) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
