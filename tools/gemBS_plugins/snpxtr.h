#ifndef SNPXTR_H_
#define SNPXTR_H_

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <htslib/vcf.h>

#include "uthash.h"

#define COMP_GZIP (1 << COMPRESS_GZIP)
#define COMP_BZIP2 (1 << COMPRESS_BZIP2)
#define COMP_XZ (1 << COMPRESS_XZ)

void error(const char *format, ...) HTS_NORETURN;

typedef struct {
  uint64_t mask;
  int n_entries;
  uint16_t *entries;
  uint8_t *name_buf;
} dbsnp_bin;

typedef struct {
  char *name;
  int min_bin;
  int max_bin;
  dbsnp_bin *bins;
  UT_hash_handle hh;
} dbsnp_ctg;

typedef struct {
  char *name;
  UT_hash_handle hh;
} snp;
  
typedef struct {
  bcf_hdr_t *hdr;
  char *output_filename;
  char *snp_filename;
  FILE *output_file;
  snp *snp_hash;
  char *dbSNP_name;
  dbsnp_ctg *dbSNP;
  char **dbSNP_prefix;
  char *dbSNP_header;
  uint16_t n_dbSNP_prefixes;
  int *gt;
  int compress;
  int pass_index;
  int gt_index;
} args_t;

#endif // SNPXTR_H_
