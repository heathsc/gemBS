/*
 * PROJECT: GEM-Tools library
 * FILE: gt_error.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_ERROR_H_
#define GT_ERROR_H_

#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdarg.h>

#include "gt_commons.h"

// Internally to Gem-tools error codes are returned as gt_status
typedef int32_t gt_status;

// Codes gt_status
#define GT_STATUS_OK 1
#define GT_STATUS_FAIL -1

// Base-name of the sources
#define GT_ERROR_BASENAME(S) \
  ({ const char* const slash=strrchr((S),'/'); \
     slash ? slash + 1 : (S); })

// Gem-tools error/output ELD-streams
extern FILE* gt_error_stream;
extern FILE* gt_log_stream;
extern FILE* gt_debug_stream;
extern bool mute_error_stream;
extern bool mute_report_stream;
// Getters/Setters ELD-streams
#define GT_DEFAULT_REPORT_STREAM stderr
GT_INLINE FILE* gt_error_get_stream();
GT_INLINE void gt_error_set_stream(FILE* const stream);
GT_INLINE FILE* gt_log_get_stream();
GT_INLINE void gt_log_set_stream(FILE* const stream);
GT_INLINE FILE* gt_debug_get_stream();
GT_INLINE void gt_debug_set_stream(FILE* const stream);
// Mute/Articulate ELD-streams
GT_INLINE void gt_mute_error_stream();
GT_INLINE void gt_mute_report_stream();
GT_INLINE void gt_articulate_error_stream();
GT_INLINE void gt_articulate_report_stream();

// Labels
#define GT_LABEL_ERROR "Error"
#define GT_LABEL_FATAL_ERROR "Fatal error"
#define GT_LABEL_DEBUG "GTDebug"
#define GT_LABEL_LOG "GTLog"
#define GT_LABEL_WARNING "Warning"

/*
 * StackTrace Printer
 */
void gt_print_stack_trace();

/*
 * Error Signal Handler
 */
void gt_handle_error_signals();

/*
 * Print ErrNo
 */
void gt_perror();

/*
 * Low-level code report functions
 *
 * E.g. EXITING ERROR (FATAL)
 *  gt_report_error_begin_block(gt_label,gt_error_name,args...) {
 *    ..Code Block..
 *  } gt_report_end_block__exit(exit_error_code,true/false);
 *
 * E.g. NOT-EXITING
 *  gt_report_error_begin_block(gt_label,gt_report_msg,args...) {
 *    ..Code Block..
 *  } gt_report_end_block();
 */
#define gt_report_error_begin_block(gt_label,gt_error_name,args...) \
  do { \
    if (!mute_error_stream) { \
      FILE* const error_stream=gt_error_get_stream(); \
      fprintf(error_stream,gt_label" (%s:%d,%s)\n "GT_ERROR_##gt_error_name"\n", \
        GT_ERROR_BASENAME(__FILE__),__LINE__,__func__, ##args); \
      fflush(error_stream); \
    }
#define gt_report_begin_block(gt_label,stream,gt_report_msg,args...) \
  do { \
    if (!mute_report_stream) { \
      FILE* const gt_stream=gt_##stream##_get_stream(); \
      fprintf(gt_stream,gt_label" (%s:%d,%s)\n "gt_report_msg"\n", \
        GT_ERROR_BASENAME(__FILE__),__LINE__,__func__, ##args); \
      fflush(gt_stream); \
    }
#define gt_report_raw_begin_block(gt_label,stream,gt_report_msg,args...) \
  do { \
    if (!mute_report_stream) { \
      FILE* const gt_stream=gt_##stream##_get_stream(); \
      fprintf(gt_stream,gt_label":: "gt_report_msg"\n",##args); \
      fflush(gt_stream); \
    }
#define gt_report__timestamp_begin_block(gt_label,stream,gt_report_msg,args...) \
  do { \
    if (!mute_report_stream) { \
      FILE* const gt_stream=gt_##stream##_get_stream(); \
      gt_tfprintf(gt_stream,gt_label":: "gt_report_msg"\n",##args); \
      fflush(gt_stream); \
    }
#define gt_report_end_block() \
  } while (0)
#define gt_report_end_block__exit(exit_code,print_stack_trace) \
    if (print_stack_trace) gt_print_stack_trace(); \
    exit(exit_code); \
  } while (0)
/*
 * GT-Error handlers
 */
#define gt_fatal_error(gt_error_name,args...) \
  gt_report_error_begin_block(GT_LABEL_FATAL_ERROR,gt_error_name,##args) \
  gt_report_end_block__exit(1,1)
#define gt_fatal_error__perror(gt_error_name,args...) \
  gt_report_error_begin_block(GT_LABEL_FATAL_ERROR,gt_error_name,##args) \
  gt_perror(); \
  gt_report_end_block__exit(1,1)
#define gt_error(gt_error_name,args...) \
  gt_report_error_begin_block(GT_LABEL_ERROR,gt_error_name,##args) \
  gt_report_end_block()
/*
 * ErrorMsg handlers
 */
#define gt_fatal_error_msg(gt_error_msg,args...) \
  gt_report_begin_block(GT_LABEL_FATAL_ERROR,error,gt_error_msg,##args) \
  gt_report_end_block__exit(1,1)
#define gt_error_msg(gt_error_msg,args...) \
  gt_report_begin_block(GT_LABEL_ERROR,error,gt_error_msg,##args) \
  gt_report_end_block()
/*
 * GT-Warning handler
 */
#define gt_warn(gt_warning_name,args...) \
  gt_report_error_begin_block(GT_LABEL_WARNING,gt_warning_name,##args) \
  gt_report_end_block()
/*
 * GT-Exception handlers (conditional error handlers)
 */
#define gt_cond_fatal_error(condition,gt_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_fatal_error(gt_error_name,##args); \
    } \
  } while (0)
#define gt_cond_fatal_error__perror(condition,gt_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_fatal_error__perror(gt_error_name,##args); \
    } \
  } while (0)
#define gt_cond_error(condition,gt_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_error(gt_error_name,##args); \
    } \
  } while (0)
/*
 * Exception handlers (conditional error handlers)
 */
#define gt_cond_fatal_error_msg(condition,error_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_fatal_error_msg(error_msg,##args); \
    } \
  } while (0)
#define gt_cond_error_msg(condition,error_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_error_msg(error_msg,##args); \
    } \
  } while (0)
/*
 * Robust checkers
 */
#ifndef GT_NO_CONSISTENCY_CHECKS
  #define gt_fatal_check(condition,gt_error_name,args...) gt_cond_fatal_error(condition,gt_error_name,##args)
  #define gt_check(condition,gt_error_name,args...) gt_cond_error(condition,gt_error_name,##args)
  #define gt_check_block(condition) if (condition)
#else
  #define gt_fatal_check(condition,gt_error_name,args...)
  #define gt_check(condition,gt_error_name,args...)
  #define gt_check_block(condition) if (0)
#endif
/*
 * Debug Utilities
 */
#ifdef GT_DEBUG
  #define gt_debug(gt_debug_msg,args...) \
    gt_report_raw_begin_block(GT_LABEL_DEBUG,debug,gt_debug_msg,##args) \
    gt_report_end_block()
  #define gt_debug_msg(gt_debug_msg,args...) \
    gt_report_begin_block(GT_LABEL_DEBUG,debug,gt_debug_msg,##args) \
    gt_report_end_block()
  #define gt_cond_debug_msg(condition,gt_debug_msg,args...) \
    do { \
      if (__builtin_expect((condition),0)){ \
        gt_debug_msg(gt_debug_msg,##args); \
      } \
    } while (0)
  #define gt_debug_msg__stack(gt_debug_msg,args...) \
    gt_debug_msg(gt_debug_msg,args...); \
    gt_print_stack_trace();
  #define gt_debug_block(condition) if (condition)
#else
  #define gt_debug_msg(gt_debug_msg,args...)
  #define gt_debug_msg__stack(gt_debug_msg,args...)
  #define gt_cond_debug_msg(condition,gt_debug_msg,args...)
  #define gt_debug_block(condition) if (0)
#endif
/*
 * Log Utilities
 */
#define gt_log(gt_log_msg,args...) \
  gt_report__timestamp_begin_block(GT_LABEL_LOG,log,gt_log_msg,##args) \
  gt_report_end_block()
#define gt_slog(gt_log_msg,args...) \
  do { \
    if (!mute_report_stream) { \
      FILE* const gt_stream=gt_log_get_stream(); \
      fprintf(gt_stream,gt_log_msg,##args); \
      fflush(gt_stream); \
    } \
  } while (0)
#define gt_cond_log(condition,gt_log_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_log(gt_log_msg,##args); \
    } \
  } while (0)

/*
 * Time Printed Formated functions
 */
gt_status gt_vtfprintf(FILE* stream,const char* format,va_list v_args);
gt_status gt_tfprintf(FILE* stream,const char* format,...);
gt_status gt_vtprintf(const char* format,va_list v_args);
gt_status gt_tprintf(const char* format,...);

/*
 * ERROR CODES/MSG
 *   #define GT_ERROR_<CODE> "<MSG>"
 */
// Library/Program errors
#define GT_ERROR_NOT_ZERO "Value Zero. Variable %s must be non-zero"
#define GT_ERROR_INVALID_VALUE "invalid Value. Variable %s must be %s"
#define GT_ERROR_POSITION_OUT_OF_RANGE "Requested position out of range"
#define GT_ERROR_POSITION_OUT_OF_RANGE_INFO "Requested position (%"PRIu64") out of range [%"PRIu64",%"PRId64"]"
#define GT_ERROR_SELECTION_NOT_IMPLEMENTED "Library error. Selection not implemented"
#define GT_ERROR_SELECTION_NOT_VALID "Library error. Selection not valid"
#define GT_ERROR_ALG_INCONSISNTENCY "Library error. Algorithmic inconsistency, check your program"
#define GT_ERROR_NOT_IMPLEMENTED "Function/Feature not implemented yet (sorry)"

// Memory errors
#define GT_ERROR_MEM_HANDLER "Could not allocate handler"
#define GT_ERROR_MEM_ALLOC "Could not allocate memory"
#define GT_ERROR_MEM_ALLOC_INFO "Could not allocate memory (%"PRIu64" requested)"
#define GT_ERROR_MEM_CALLOC_INFO "Could not allocate memory (%"PRIu64"x%"PRIu64" requested)"
#define GT_ERROR_MEM_ALLOC_DISK "Requested %"PRIu64" Bytes. Resorting to disk"
#define GT_ERROR_MEM_ALLOC_MMAP_DISK_FAIL "Requested %"PRIu64" Bytes. Failed to mmap memory to '%s'"
#define GT_ERROR_MEM_ALLOC_MMAP_FAIL "Requested %"PRIu64" Bytes. Failed to mmap memory"
#define GT_ERROR_MEM_REALLOC "Could not reallocate memory"
#define GT_ERROR_MEM_CURSOR_OUT_OF_SEGMENT "Current memory cursor is out of boundaries (Segmentation fault)"
#define GT_ERROR_MEM_CURSOR_SEEK "Could not seek to address %"PRIu64". Out of boundaries (Segmentation fault)"
#define GT_ERROR_MEM_ALG_FAILED "Failed aligning the memory address to the specified boundary"
#define GT_ERROR_NULL_HANDLER "Null handler or fields not properly allocated"
#define GT_ERROR_NULL_HANDLER_INFO "Null handler %s "

// System errors
#define GT_ERROR_SYS_ERROR "System error signal raised"
#define GT_ERROR_SYS_MMAP "Mmap call error"
#define GT_ERROR_SYS_MMAP_FILE "Could not map file '%s' to memory"
#define GT_ERROR_SYS_UNMAP "Could not unmap memory"
#define GT_ERROR_SYS_THREAD "Could not create thread"
#define GT_ERROR_SYS_PIPE "Could not create pipe"
#define GT_ERROR_SYS_MUTEX "Mutex call error"
#define GT_ERROR_SYS_MUTEX_INIT "Mutex initialization error"
#define GT_ERROR_SYS_MUTEX_DESTROY "Mutex destroy call error"
#define GT_ERROR_SYS_COND_VAR "Conditional variable call error"
#define GT_ERROR_SYS_COND_VAR_INIT "Conditional variable initialization error"
#define GT_ERROR_SYS_COND_VAR_DESTROY "Conditional variable destroy call error"
#define GT_ERROR_SYS_MKSTEMP "Could not create a temporal file (mkstemp:'%s')"
#define GT_ERROR_SYS_HANDLE_TMP "Failed to handle temporal file"
#define GT_ERROR_SYS_OPEN_PIPE "Failed to open pipe"
#define GT_ERROR_SYS_CLOSE_PIPE "Failed to close pipe"
#define GT_ERROR_SYS_FORK "Failed to fork"
#define GT_ERROR_SYS_EXEC "Failed to exec %s %s"

// String errors
#define GT_ERROR_STRING_STATIC "Could not perform operation on static string"

// File errors
#define GT_ERROR_FILE_STAT "Could not stat file '%s'"
#define GT_ERROR_FILE_OPEN "Could not open file '%s'"
#define GT_ERROR_FILE_READ "Could not read from file '%s'"
#define GT_ERROR_FILE_READ_FD "Could not read from file descriptor"
#define GT_ERROR_FILE_WRITE "Could not write to file '%s'"
#define GT_ERROR_FILE_SEEK "Could not seek file '%s' to %"PRIu64" "
#define GT_ERROR_FILE_CLOSE "Could not close file '%s'"
#define GT_ERROR_FILE_FORMAT "Could not determine file format"
#define GT_ERROR_FILE_GZIP_OPEN "Could not open GZIPPED file '%s'"
#define GT_ERROR_FILE_GZIP_NO_ZLIB "Could not open GZIPPED file '%s': no zlib support compiled in"
#define GT_ERROR_FILE_BZIP2_OPEN "Could not open BZIPPED file '%s'"
#define GT_ERROR_FILE_BZIP2_NO_BZLIB "Could not open BZIPPED file '%s': no bzlib support compiled in"
#define GT_ERROR_FILE_FDOPEN "Could not fdopen file descriptor"
#define GT_ERROR_FILE_XZIP_OPEN "Could not open XZIPPED file '%s'"

// Pipe errors
#define GT_ERROR_PIPE_BAD_PIPE "Invalid pipe command: '%s'"
#define GT_ERROR_PIPE_BIDIRECTIONAL_PIPE "Can not use bidirectional pipes: '%s'"
#define GT_ERROR_PIPE_NOT_READ "Pipe is not for reading: '%s'"
#define GT_ERROR_PIPE_NOT_WRITE "Pipe is not for writing: '%s'"

// Output errors
#define GT_ERROR_FPRINTF "Printing output. 'fprintf' call failed"
#define GT_ERROR_SPRINTF "Printing output. 'sprintf' call failed"
#define GT_ERROR_TPRINTF "Printing output. 'tprintf' call failed"
#define GT_ERROR_BPRINTF "Printing output. Buffer print formated 'gt_bprintf' call failed"
#define GT_ERROR_OFPRINTF "Printing output. Output File print formated 'gt_ofprintf' call failed"
#define GT_ERROR_BOFPRINTF "Printing output. Buffered Output file print formated 'gt_bofprintf' call failed"

#define GT_ERROR_PRINT_FORMAT "Incorrect print format. Expected format character"

// Template/Alignment/Map/Misms errors
#define GT_ERROR_MISMS_TYPE "Misms incorrect type"
#define GT_ERROR_MISMS_TRANSITION "Incorrect mismatch transition (Same base at both sides)"
#define GT_ERROR_MISMS_SPLICE_POS "Splicing distance must be positive (non-zero)"
#define GT_ERROR_COUNTERS_POS_STRATUM "Stratum must be strictly positive (stratum>0)"
#define GT_ERROR_MAP_MISMS_NOT_PARSED "Map's mismatches not parsed yet"
#define GT_ERROR_MAP_NEG_LENGTH "Negative Map total length"
#define GT_ERROR_MAP_NEG_MAPPED_BASES "Negative number of bases mapped"
#define GT_ERROR_ALIGNMENT_READ_QUAL_LENGTH "Read and quality length differs"
#define GT_ERROR_ALIGNMENT_MAPS_NOT_PARSED "Alignment's maps not parsed yet"
#define GT_ERROR_ALIGNMENT_INCONSISTENT_COUNTERS "Alignment inconsistency. Maps inconsistent with counters values"
#define GT_ERROR_ALIGNMENT_NOT_SCORED "Alignments have no valid score.  MAPQ can not be calculated unless the map file has been processed using score_reads"
#define GT_ERROR_TEMPLATE_ZERO_BLOCKS "Zero alignment blocks (num_blocks_template==0)"
#define GT_ERROR_TEMPLATE_TOO_MANY_BLOCKS "Template contains already two ends"
#define GT_ERROR_TEMPLATE_BLOCKS_EXCESS "Template blocks exceeds two ends"
#define GT_ERROR_TEMPLATE_MMAP_NULL "Template mmap is null (all maps are null)"
#define GT_ERROR_TEMPLATE_INCONSISTENT_MMAPS_ALIGNMENT "Template inconsistency. Multimaps' members must be contained by single alignments"
#define GT_ERROR_TEMPLATE_INCONSISTENT_MMAPS_ATTRB_RELATION "Template inconsistency. Incorrect number of mmaps and mmaps' attributes"
#define GT_ERROR_TEMPLATE_INCONSISTENT_NUM_MAPS_RELATION "Template inconsistency. Incorrect number of matches' elements (check num_blocks_template)"
#define GT_ERROR_TEMPLATE_INCONSISTENT_NUM_BLOCKS "Template inconsistency. Number of blocks must be the same across templates"
#define GT_ERROR_TEMPLATE_INCONSISTENT_COUNTERS "Template inconsistency. MMaps inconsistent with counters values"
#define GT_ERROR_TEMPLATE_ADD_BAD_NUM_BLOCKS "Trying to add wrong number of blocks to the template"
#define GT_ERROR_PALIGN_BAD_NUM_BLOCKS "Invalid Paired-alignment. Wrong number of alignment blocks (%"PRIu64")"
#define GT_ERROR_TEMPLATE_NOT_SCORED "Alignments have no valid score.  MAPQ can not be calculated unless the map file has been processed using score_reads"
	  
// Sequence Archive/Segmented Sequence errors
#define GT_ERROR_SEGMENTED_SEQ_IDX_OUT_OF_RANGE "Error accessing segmented sequence. Index %"PRIu64" out out range [0,%"PRIu64")"
#define GT_ERROR_CDNA_IT_OUT_OF_RANGE "Error seeking sequence. Index %"PRIu64" out out range [0,%"PRIu64")"
#define GT_ERROR_SEQ_ARCHIVE_WRONG_TYPE "Wrong sequence archive type"
#define GT_ERROR_SEQ_ARCHIVE_NOT_FOUND "Sequence '%s' not found in reference archive"
#define GT_ERROR_SEQ_ARCHIVE_POS_OUT_OF_RANGE "Requested position '%"PRIu64"' out of sequence boundaries"
#define GT_ERROR_SEQ_ARCHIVE_CHUNK_OUT_OF_RANGE "Requested sequence string [%"PRIu64",%"PRIu64") out of sequence '%s' boundaries"
#define GT_ERROR_GEMIDX_SEQ_ARCHIVE_NOT_FOUND "GEMIdx. Sequence '%s' not found in reference archive"
#define GT_ERROR_GEMIDX_INTERVAL_NOT_FOUND "GEMIdx. Interval relative to sequence '%s' not found in reference archive"

// Stats vector
#define GT_ERROR_VSTATS_INVALID_MIN_MAX "Invalid step range for stats vector, min_value <= max_value"

/*
 * Parsing FASTQ File format errors
 */
// IFP (Input FASTA Parser). General
#define GT_ERROR_PARSE_FASTA "Parsing FASTA/FASTQ error(%s:%"PRIu64":%"PRIu64")"

/*
 * Parsing MAP File format errors
 */
// IMP (Input MAP Parser). General
#define GT_ERROR_PARSE_MAP "Parsing MAP error(%s:%"PRIu64")"
#define GT_ERROR_PARSE_MAP_BAD_FILE_FORMAT "Parsing MAP error(%s:%"PRIu64"). Not a MAP file"
#define GT_ERROR_PARSE_MAP_BAD_NUMBER_FIELDS "Parsing MAP error(%s:%"PRIu64"). Wrong number of TAB separated fields (%"PRIu64")"
#define GT_ERROR_PARSE_MAP_BAD_READ_QUAL_LENGTH "Parsing MAP error(%s:%"PRIu64"). Mismatching Read length (%"PRIu64") and Quality length (%"PRIu64")"
#define GT_ERROR_PARSE_MAP_COUNTERS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Error parsing counters"
#define GT_ERROR_PARSE_MAP_BAD_TEMPLATE_SEP "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Read character '%c' not valid (%s)"
#define GT_ERROR_PARSE_MAP_DIFF_TEMPLATE_BLOCKS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Different number of template blocks {read(%"PRIu64"),qualities(%"PRIu64")}"
#define GT_ERROR_PARSE_MAP_NOT_AN_ALIGNMENT "Parsing MAP error(%s:%"PRIu64"). File doesn't contains simple alignments (use template)"
#define GT_ERROR_PARSE_MAP_MISMS_ALREADY_PARSED "Parsing MAP error(%s:%"PRIu64"). Mismatch string already parsed or null lazy-parsing handler"
#define GT_ERROR_PARSE_MAP_NOT_IMPLEMENTED "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Feature not implemented yet (sorry)"
#define GT_ERROR_PARSE_MAP_PREMATURE_EOL "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Premature End-of-line found"
// IMP (Input MAP Parser). Parsing Read Errors
#define GT_ERROR_PARSE_MAP_READ_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing read, bad character found"
// IMP (Input MAP Parser). Parsing Qualities Errors
#define GT_ERROR_PARSE_MAP_QUAL_BAD_SEPARATOR "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing quality string, bad block-separator found"
#define GT_ERROR_PARSE_MAP_QUAL_BAD_LENGTH "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing quality string, wrong length (w.r.t. read length)"
#define GT_ERROR_PARSE_MAP_QUAL_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing quality string, wrong quality value (bad character)"
// IMP (Input MAP Parser). Parsing Counters Errors
#define GT_ERROR_PARSE_MAP_COUNTERS_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing counters, bad character found"
// IMP (Input MAP Parser). Parsing Maps Errors
#define GT_ERROR_PARSE_MAP_BAD_NUMBER_OF_BLOCKS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing maps, wrong number of blocks"
#define GT_ERROR_PARSE_MAP_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing maps, bad character found"
#define GT_ERROR_PARSE_MAP_INCONSISTENT_BLOCKS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing maps, block(s) length doesn't match the read length"
#define GT_ERROR_PARSE_MAP_SPLIT_MAP_BAD_NUM_ACCEPTORS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing split-map, bad number of acceptors"
#define GT_ERROR_PARSE_MAP_SPLIT_MAP_BAD_NUM_DONORS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing split-map, bad number of donors"
// IMP (Input MAP Parser). Parsing Mismatch String Errors
#define GT_ERROR_PARSE_MAP_MISMS_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing mismatch string, bad character found"
#define GT_ERROR_PARSE_MAP_MISMS_BAD_MISMS_POS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing mismatch string, unsorted mismatches"

/*
 * Parsing SAM File format errors
 */
// ISP (Input SAM Parser). General
#define GT_ERROR_PARSE_SAM "Parsing SAM error(%s:%"PRIu64":%"PRIu64")"
#define GT_ERROR_PARSE_SAM_HEADER_NOT_SAM "Parsing SAM Header error(%s). Header file not in SAM format"
#define GT_ERROR_PARSE_SAM_HEADER_MISSING_TAG "Parsing SAM Header error. @%s record missing %s tag"
#define GT_ERROR_PARSE_SAM_HEADER_DUPLICATE_TAG "Parsing SAM Header error. @%s record with duplicate %s tag (%s)"
#define GT_ERROR_PARSE_SAM_BAD_FILE_FORMAT "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Not a SAM file"
#define GT_ERROR_PARSE_SAM_BAD_CHARACTER "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Bad character found"
#define GT_ERROR_PARSE_SAM_UNMAPPED_XA "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Unmapped read contains XA field (inconsistency)"
#define GT_ERROR_PARSE_SAM_PREMATURE_EOL "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Premature End-of-line found"
#define GT_ERROR_PARSE_SAM_EXPECTED_NUMBER "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Expected number"
#define GT_ERROR_PARSE_SAM_WRONG_READ_CONTENT "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Read in template doesn't match previously parse reads with same tag"
#define GT_ERROR_PARSE_SAM_CIGAR_PREMATURE_END "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Premature end of CIGAR string"
#define GT_ERROR_PARSE_SAM_WRONG_NUM_XA "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Wrong number of eXtra mAps (as to pair them)"
#define GT_ERROR_PARSE_SAM_UNSOLVED_PENDING_MAPS "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Failed to pair maps"
#define GT_ERROR_PARSE_SAM_DUPLICATE_SEQUENCE_TAG "Parsing SAM error.  Sequence tag "PRIgts" not unique when pairing unique mapping reads"

#define GT_ERROR_SAM_OUTPUT_UNKNOWN_RG_ID "SAM output error.  Read group ID %s not found in SAM headers"
#define GT_ERROR_SAM_OUTPUT_NO_HEADER_FOR_RG "SAM output error.  No SAM header was specified, so read group ID can not be matched"

/*
 * Output File
 */
#define GT_ERROR_OUTPUT_FILE_INCONSISTENCY "Output file state inconsistent"
#define GT_ERROR_OUTPUT_FILE_FAIL_WRITE "Output file. Error writing to to file"
#define GT_ERROR_BUFFER_SAFETY_DUMP "Output buffer. Could not perform safety dump"

#define GT_ERROR_OUTPUT_SAM_NO_PRIMARY_ALG "Output SAM. No primary alignment specified"

/*
 * Map Alignment
 */
#define GT_ERROR_MAP_ALG_WRONG_ALG "(Re)Aligning Map. Wrong alignment"
#define GT_ERROR_MAP_RECOVER_MISMS_WRONG_BASE_ALG "Recovering mismatches from map. Wrong initial alignment"

/*
 * General purpose checkers
 */
//#define GT_NULL_CHECK(object) gt_fatal_check((object)==NULL,NULL_HANDLER_INFO,((char*)GT_QUOTE(object)))
#define GT_ZERO_CHECK(object) gt_fatal_check((object)==0,NOT_ZERO,((char*)GT_QUOTE(object)))
#define GT_INVALID_CASE() gt_fatal_error(SELECTION_NOT_VALID)
/* Eclipse debugging definitions */
#define GT_NULL_CHECK(object) if ((object)==NULL) { printf("%d\n",(*(int*)object)); }

#endif /* GT_ERROR_H_ */
