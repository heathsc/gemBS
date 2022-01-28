/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_file.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#define _GNU_SOURCE

#include <sys/wait.h>
#ifdef ZLIB
#include <zlib.h>
#endif
#ifdef ZLIB
#include <bzlib.h>
#endif
#include "gt_pipe_io.h"
#include "gt_input_file.h"

// Internal constants
#define GT_INPUT_BUFFER_SIZE GT_BUFFER_SIZE_64M

/*
 * Basic I/O functions
 */

gt_input_file* gt_input_stream_general_open(FILE* stream,gt_file_format format) {
  GT_NULL_CHECK(stream);
  // Allocate handler
  gt_input_file* input_file = gt_alloc(gt_input_file);
  // Input file
  input_file->file_name = GT_STREAM_FILE_NAME;
  input_file->file_type = STREAM;
  input_file->file = stream;
  input_file->fildes = -1;
  input_file->eof = feof(stream);
  input_file->file_size = UINT64_MAX;
  gt_cond_fatal_error(pthread_mutex_init(&input_file->input_mutex, NULL),SYS_MUTEX_INIT);
  // Auxiliary Buffer (for synch purposes)
  input_file->file_buffer = gt_malloc(GT_INPUT_BUFFER_SIZE);
  input_file->buffer_size = 0;
  input_file->buffer_begin = 0;
  input_file->buffer_pos = 0;
  input_file->global_pos = 0;
  input_file->processed_lines = 0;
  // ID generator
  input_file->processed_id = 0;
  return input_file;
}

gt_input_file* gt_input_file_pipe_open(char *const file_name,gt_file_format format)
{
  char *trimmed_name=NULL;
  gt_cond_fatal_error(gt_pipe_io_check_command(file_name,&trimmed_name)!=GT_PIPE_IO_READ,PIPE_NOT_READ,file_name);
  GT_NULL_CHECK(trimmed_name);
  int fd=gt_pipe_io_child_open(GT_PIPE_IO_READ,trimmed_name);
  gt_cond_fatal_error(fd<0,SYS_PIPE); // Should never happen
  FILE *stream=fdopen(fd,"r");
  gt_cond_fatal_error(stream==NULL,FILE_OPEN,trimmed_name); // Should never happen
  free(trimmed_name);
  return gt_input_stream_general_open(stream,format);
}

gt_input_file* gt_input_file_general_open(char* const file_name,const bool mmap_file,gt_file_format format)
{
  GT_NULL_CHECK(file_name);
  // Check for pipe
  if(strchr(file_name,'|')) return gt_input_file_pipe_open(file_name,format);
  // Allocate handler
  gt_input_file* input_file = gt_alloc(gt_input_file);
  // Input file
  struct stat stat_info;
  unsigned char tbuf[6];
  int i;
  gt_cond_fatal_error(stat(file_name,&stat_info)==-1,FILE_STAT,file_name);
  input_file->file_name = file_name;
  input_file->file_size = stat_info.st_size;
  input_file->eof = (input_file->file_size==0);
  gt_cond_fatal_error(pthread_mutex_init(&input_file->input_mutex,NULL),SYS_MUTEX_INIT);
  if (mmap_file) {
    input_file->file = NULL;
    input_file->fildes = open(file_name,O_RDONLY,0); // TODO: O_NOATIME condCompl (Thanks Jordi Camps)
    gt_cond_fatal_error(input_file->fildes==-1,FILE_OPEN,file_name);
    input_file->file_buffer =
      (uint8_t*) mmap(0,input_file->file_size,PROT_READ,MAP_PRIVATE,input_file->fildes,0);
    gt_cond_fatal_error(input_file->file_buffer==MAP_FAILED,SYS_MMAP_FILE,file_name);
    input_file->file_type = MAPPED_FILE;
  } else {
    input_file->fildes = -1;
    gt_cond_fatal_error(!(input_file->file=fopen(file_name,"r")),FILE_OPEN,file_name);
    input_file->file_type = REGULAR_FILE;
    if(S_ISREG(stat_info.st_mode)) {
      // Regular file - check if gzip or bzip compressed
      i=(int)fread(tbuf,(size_t)1,(size_t)6,input_file->file);
		if(i == 6) {
			if(tbuf[0]==0x1f && tbuf[1]==0x8b && tbuf[2]==0x08) {
				input_file->file_type=GZIPPED_FILE;
				fclose(input_file->file);
#ifdef HAVE_ZLIB
				gt_cond_fatal_error(!(input_file->file=(void *)gzopen(file_name,"r")),FILE_GZIP_OPEN,file_name);
#else
				gt_fatal_error(FILE_GZIP_NO_ZLIB,file_name);
#endif
			} else if(tbuf[0]=='B' && tbuf[1]=='Z' && tbuf[2]=='h' && tbuf[3]>='0' && tbuf[3]<='9') {
				fseek(input_file->file,0L,SEEK_SET);
				input_file->file_type=BZIPPED_FILE;
#ifdef HAVE_BZLIB
				input_file->file=BZ2_bzReadOpen(&i,input_file->file,0,0,NULL,0);
				gt_cond_fatal_error(i!=BZ_OK,FILE_BZIP2_OPEN,file_name);
#else
				gt_fatal_error(FILE_BZIP2_NO_BZLIB,file_name);
#endif
			} else  if (tbuf[0] == 0xfd && tbuf[1] == '7' && tbuf[2] == 'z' && tbuf[3] == 'X' && tbuf[4] == 'Z' && tbuf[5] == 0) {
				fseek(input_file->file,0L,SEEK_SET);
				input_file->file_type=XZIPPED_FILE;
				char *nbuf = NULL;
				asprintf(&nbuf, "xzcat %s", file_name);
				gt_cond_fatal_error(nbuf == NULL, MEM_ALLOC);
				input_file->file=popen(nbuf, "r");
//				fprintf(stderr, "xzip:  cmd line = '%s', file = %p\n", nbuf, input_file->file);
				gt_cond_fatal_error(input_file == NULL, FILE_XZIP_OPEN, file_name);
				free(nbuf);
			} else {
				fseek(input_file->file,0L,SEEK_SET);
			}
      } else {
        fseek(input_file->file,0L,SEEK_SET);
      }
    } else {
      input_file->eof=0;
    }
    input_file->file_buffer = gt_malloc(GT_INPUT_BUFFER_SIZE);
  }
  // Auxiliary Buffer (for synch purposes)
  input_file->buffer_size = 0;
  input_file->buffer_begin = 0;
  input_file->buffer_pos = 0;
  input_file->global_pos = 0;
  input_file->processed_lines = 0;
  // ID generator
  input_file->processed_id = 0;
  return input_file;
}
/*
 * POST: Closes the gt_input_file
 * RETURN VALUE: Returns zero on success and error code
 */
gt_status gt_input_file_close(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  gt_status status = GT_INPUT_FILE_OK;
#ifdef HAVE_BZLIB
  int bzerr;
#endif
  switch (input_file->file_type) {
    case REGULAR_FILE:
      gt_free(input_file->file_buffer);
      if (fclose(input_file->file)) status = GT_INPUT_FILE_CLOSE_ERR;
      break;
    case GZIPPED_FILE:
      gt_free(input_file->file_buffer);
#ifdef HAVE_ZLIB
      if (gzclose((gzFile)input_file->file)) status = GT_INPUT_FILE_CLOSE_ERR;
#endif
      break;
    case BZIPPED_FILE:
      gt_free(input_file->file_buffer);
#ifdef HAVE_BZLIB
      BZ2_bzReadClose(&bzerr,input_file->file);
      if (bzerr!=BZ_OK) status = GT_INPUT_FILE_CLOSE_ERR;
#endif
      break;
	case XZIPPED_FILE:
	  gt_free(input_file->file_buffer);
	  pclose(input_file->file);
	  break;
    case MAPPED_FILE:
      gt_cond_error(munmap(input_file->file,input_file->file_size)==-1,SYS_UNMAP);
      if (close(input_file->fildes)) status = GT_INPUT_FILE_CLOSE_ERR;
      break;
    case STREAM:
      gt_free(input_file->file_buffer);
      // In case the stream is from a pipe where we forked a shell we should wait for the child to die
      signal(SIGCHLD,SIG_DFL);
      int i;
      while(waitpid(-1,&i,WNOHANG)>0);
      break;
  }
  gt_free(input_file);
  return status;
}

/*
 * Accessors (Mutex,ID,...) functions
 */
GT_INLINE void gt_input_file_lock(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  gt_cond_fatal_error(pthread_mutex_lock(&input_file->input_mutex),SYS_MUTEX);
}
GT_INLINE void gt_input_file_unlock(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  gt_cond_fatal_error(pthread_mutex_unlock(&input_file->input_mutex),SYS_MUTEX);
}
GT_INLINE uint64_t gt_input_file_next_id(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  return (input_file->processed_id)++;
}

/*
 * Basic line functions
 */
GT_INLINE size_t gt_input_file_dump_to_buffer(gt_input_file* const input_file,gt_vector* const buffer_dst) { // FIXME: If mmap file, internal buffer is just pointers to mem
  GT_INPUT_FILE_CHECK(input_file);
  // Copy internal file buffer to buffer_dst
  const uint64_t chunk_size = input_file->buffer_pos-input_file->buffer_begin;
  if (gt_expect_false(chunk_size==0)) return 0;
  gt_vector_reserve_additional(buffer_dst,chunk_size);
  memcpy(gt_vector_get_mem(buffer_dst,uint8_t)+gt_vector_get_used(buffer_dst),
      input_file->file_buffer+input_file->buffer_begin,chunk_size);
  gt_vector_add_used(buffer_dst,chunk_size);
  // Update position
  input_file->buffer_begin=input_file->buffer_pos;
  // Return number of written bytes
  return chunk_size;
}
GT_INLINE size_t gt_input_file_fill_buffer(gt_input_file* const input_file) {
#ifdef HAVE_BZLIB
  int bzerr;
#endif
  GT_INPUT_FILE_CHECK(input_file);
  input_file->global_pos += input_file->buffer_size;
  input_file->buffer_pos = 0;
  input_file->buffer_begin = 0;
  if (gt_expect_true(
      (input_file->file_type==STREAM && !feof(input_file->file)) ||
      (input_file->file_type==XZIPPED_FILE && !feof(input_file->file)) ||
      (input_file->file_type==REGULAR_FILE && !feof(input_file->file)))) {
    input_file->buffer_size =
        fread(input_file->file_buffer,sizeof(uint8_t),GT_INPUT_BUFFER_SIZE,input_file->file);
    if (input_file->buffer_size==0) {
      input_file->eof = true;
    }
    return input_file->buffer_size;
  } else if (input_file->file_type==MAPPED_FILE && input_file->global_pos < input_file->file_size) {
    input_file->buffer_size = input_file->file_size-input_file->global_pos;
    return input_file->buffer_size;
#ifdef HAVE_ZLIB
  } else if (input_file->file_type==GZIPPED_FILE && !gzeof((gzFile)input_file->file)) {
    input_file->buffer_size = gzread((gzFile)input_file->file,input_file->file_buffer,GT_INPUT_BUFFER_SIZE);
    if (input_file->buffer_size==0) {
      input_file->eof = true;
    }
    return input_file->buffer_size;
#endif
#ifdef HAVE_BZLIB
  } else if (input_file->file_type==BZIPPED_FILE) {
    input_file->buffer_size = BZ2_bzRead(&bzerr,input_file->file,input_file->file_buffer,GT_INPUT_BUFFER_SIZE);
    if(input_file->buffer_size==0) {
      input_file->eof=true;
    }
    return input_file->buffer_size;
#endif
  } else {
    input_file->eof = true;
    return 0;
  }
}
GT_INLINE size_t gt_input_file_next_line(gt_input_file* const input_file,gt_vector* const buffer_dst) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  GT_INPUT_FILE_CHECK_BUFFER__DUMP(input_file,buffer_dst);
  if (input_file->eof) return GT_INPUT_FILE_EOF;
  // Read line
  while (gt_expect_true(!input_file->eof &&
      GT_INPUT_FILE_CURRENT_CHAR(input_file)!=EOL &&
      GT_INPUT_FILE_CURRENT_CHAR(input_file)!=DOS_EOL)) {
    GT_INPUT_FILE_NEXT_CHAR__DUMP(input_file,buffer_dst);
  }
  // Handle EOL
  GT_INPUT_FILE_HANDLE_EOL(input_file,buffer_dst);
  return GT_INPUT_FILE_LINE_READ;
}
GT_INLINE size_t gt_input_file_next_record(
    gt_input_file* const input_file,gt_vector* const buffer_dst,gt_string* const first_field,
    uint64_t* const num_blocks,uint64_t* const num_tabs) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  GT_INPUT_FILE_CHECK_BUFFER__DUMP(input_file,buffer_dst);
  if (input_file->eof) return GT_INPUT_FILE_EOF;
  // Read line
  uint64_t const begin_line_pos_at_file = input_file->buffer_pos;
  uint64_t const begin_line_pos_at_buffer = gt_vector_get_used(buffer_dst);
  uint64_t current_pfield = 0, length_first_field = 0;
  while (gt_expect_true(!input_file->eof &&
      GT_INPUT_FILE_CURRENT_CHAR(input_file)!=EOL &&
      GT_INPUT_FILE_CURRENT_CHAR(input_file)!=DOS_EOL)) {
    if (current_pfield==0) {
      if (gt_expect_false(GT_INPUT_FILE_CURRENT_CHAR(input_file)==TAB)) {
        ++current_pfield; ++(*num_tabs);
      } else {
        ++length_first_field;
      }
    } else if (current_pfield==1) {
      if (gt_expect_false(GT_INPUT_FILE_CURRENT_CHAR(input_file)==SPACE)) {
        ++(*num_blocks);
      } else if (gt_expect_false(GT_INPUT_FILE_CURRENT_CHAR(input_file)==TAB)) {
        ++current_pfield; ++(*num_tabs); ++(*num_blocks);
      }
    } else {
      if (gt_expect_false(GT_INPUT_FILE_CURRENT_CHAR(input_file)==TAB)) {
        ++current_pfield; ++(*num_tabs);
      }
    }
    GT_INPUT_FILE_NEXT_CHAR__DUMP(input_file,buffer_dst);
  }
  // Handle EOL
  GT_INPUT_FILE_HANDLE_EOL(input_file,buffer_dst);
  // Set first field (from the input_file_buffer or the buffer_dst)
  if (first_field) {
    char* first_field_begin;
    if (input_file->buffer_pos <= begin_line_pos_at_file) {
      gt_input_file_dump_to_buffer(input_file,buffer_dst); // Forced to dump to buffer
      first_field_begin = gt_vector_get_elm(buffer_dst,begin_line_pos_at_buffer,char);
    } else {
      first_field_begin = (char*)input_file->file_buffer + begin_line_pos_at_file;
    }
    gt_string_set_nstring(first_field,first_field_begin,length_first_field);
  }
  return GT_INPUT_FILE_LINE_READ;
}

#define GT_INPUT_FILE_TEST_NEXT_CHAR(input_file,buffer_centinel) \
  ++buffer_centinel; \
  if (gt_expect_false(buffer_centinel >= input_file->buffer_size)) return true;
#define GT_INPUT_FILE_TEST_CURRENT_CHAR(input_file,buffer_centinel) input_file->file_buffer[buffer_centinel]
GT_INLINE bool gt_input_file_next_record_cmp_first_field(gt_input_file* const input_file,gt_string* const first_field) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_STRING_CHECK(first_field);
  if (gt_expect_false(input_file->eof || input_file->buffer_pos >= input_file->buffer_size)) return true;
  // Read line
  char* const tag_begin = (char*)(input_file->file_buffer+input_file->buffer_pos);
  uint64_t buffer_centinel = input_file->buffer_pos;
  while (gt_expect_true(!input_file->eof &&
      GT_INPUT_FILE_TEST_CURRENT_CHAR(input_file,buffer_centinel)!=EOL &&
      GT_INPUT_FILE_TEST_CURRENT_CHAR(input_file,buffer_centinel)!=DOS_EOL)) {
    if (gt_expect_false(
        (GT_INPUT_FILE_TEST_CURRENT_CHAR(input_file,buffer_centinel)==SPACE ||
         GT_INPUT_FILE_TEST_CURRENT_CHAR(input_file,buffer_centinel)==TAB) )) {
      char* const tag_end = (char*)(input_file->file_buffer+buffer_centinel);
      uint64_t tag_lenth = tag_end-tag_begin;
      if (tag_lenth>2 && tag_begin[tag_lenth-2]==SLASH) tag_lenth-=2;
      if (first_field->length != tag_lenth) return false;
      return gt_strneq(first_field->buffer,tag_begin,tag_lenth);
    }
    GT_INPUT_FILE_TEST_NEXT_CHAR(input_file,buffer_centinel);
  }
  return true;
}

/*
 * Line Readers (thread-unsafe, must call mutex functions before)
 */
GT_INLINE uint64_t gt_input_file_add_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  // Read lines
  uint64_t lines_read = 0;
  while (lines_read<num_lines && gt_input_file_next_line(input_file,buffer_dst)) {
    ++lines_read;
  }
  // Dump remaining content into the buffer
  gt_input_file_dump_to_buffer(input_file,buffer_dst);
  if (lines_read > 0 && *gt_vector_get_last_elm(buffer_dst,char) != EOL) {
    gt_vector_insert(buffer_dst,EOL,char);
  }
  input_file->processed_lines+=lines_read;
  return lines_read;
}
GT_INLINE uint64_t gt_input_file_get_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  // Clear dst buffer
  gt_vector_clear(buffer_dst);
  // Read lines
  return gt_input_file_add_lines(input_file,buffer_dst,num_lines);
}
