/*
 * dbSNP_struct.h
 *
 *  Created on: Jan 23, 2020
 *      Author: heath
 */
#include <stdlib.h>
#include <inttypes.h>
#include "uthash.h"

#ifndef INCLUDE_DBSNP_H_
#define INCLUDE_DBSNP_H_

typedef struct {
	uint64_t mask;
	uint64_t fq_mask;
	int n_entries;
	uint16_t *entries;
	uint8_t *name_buf;
} dbsnp_bin_t;

typedef struct {
	char *name;
	int min_bin;
	int max_bin;
	uint64_t file_offset;
	dbsnp_bin_t *bins;
	UT_hash_handle hh;
} dbsnp_ctg_t;

typedef struct {
	char *filename;
	FILE *fp;
	dbsnp_ctg_t *dbSNP;
	uint16_t n_dbSNP_prefixes;
	size_t dbSNP_bufsize;
	char **dbSNP_prefix;
	char *dbSNP_header;
} dbsnp_header_t;

dbsnp_header_t * load_dbSNP_header(char * const filename);
bool load_dbSNP_ctg(const dbsnp_header_t * const hdr, dbsnp_ctg_t * const ctg);
void unload_dbSNP_ctg(dbsnp_ctg_t * const ctg);
uint8_t dbSNP_lookup_name(const dbsnp_header_t *const hdr, const dbsnp_ctg_t * ctg, char * const rs, size_t * const rs_len, const uint32_t x);

#endif /* INCLUDE_DBSNP_H_ */
