//
//  dbSNP_idx.h
//  dbSNP
//
//  Created by Simon Heath on 15/11/2017.
//  Copyright 2017 Simon Heath. All rights reserved.
//

#ifndef dbSNP_idx_h
#define dbSNP_idx_h

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <assert.h>

#include "mm.h"
#include "uthash.h"


#define INITIAL_ENTRY_SIZE 4
#define INITIAL_NAME_BUF_SIZE 32
#define ITEMS_PER_BLOCK 2048
#define IDX_MAGIC 0xd7278434

typedef struct {
	void *mem;
	size_t size;
	size_t len;
} buffer_t;

typedef struct {
   uint64_t mask;
   uint64_t fq_mask;
   uint16_t name_buf_size;
   uint16_t name_buf_idx;
   uint8_t n_entries;
   uint8_t entry_size;
   uint16_t *entry;
   unsigned char *name_buf;
} bin;

typedef struct {
	char *name;
	FILE *fp;
	bool sorted;
	bool filter;
	bool processed;
	bool read;
	UT_hash_handle hh;
} file_t;

typedef struct _contig {
	struct _contig *next;
	char *name;
	char *rname;
	int name_len;
	bin *bins;
	uint32_t min_bin;
	uint32_t max_bin;
	file_t *first_file;
	uint64_t offset;
	bool in_queue;
	pthread_mutex_t mut;
	UT_hash_handle hh;
} contig;

typedef struct {
	uint16_t ix;
	char *pref;
	UT_hash_handle hh;
} prefix;

typedef struct {
	char *name;
	char *alias;
	UT_hash_handle hh;
} chrom_alias_t;

typedef enum { dbsnp_auto, dbsnp_bed, dbsnp_json, dbsnp_vcf } dbsnp_input_type_t;
typedef enum { block_empty, block_uncompressed, block_active, block_compressed } block_state_t;

typedef struct {
	buffer_t cblock;
	buffer_t ublock;
	uint64_t *off_ptr;
	block_state_t state;
	uint8_t flag;
} comp_block_t;

typedef struct {
	char *name;
	char *cname;
	uint32_t pos;
	int name_len;
	int name_off;
	int cname_len;
	double maf;
	bool ok;
} snp_t;

typedef struct {
	char *output_file;
	char *chrom_alias_file;
	char *select_file;
	void *select_hash;
	FILE *outfile;
	int threads;
	int read_jobs;
	int compress_jobs;
	dbsnp_input_type_t input_type;
	char *header;
	contig *contigs;
	contig *contig_queue;
	chrom_alias_t *aliases;
	prefix *prefixes;
	file_t *sorted;
	file_t *unsorted;
	file_t *used_files;
	buffer_t *output_buf;
	uint64_t n_snps;
	uint64_t n_snps_maf_filtered;
	uint64_t max_buf_size;
	comp_block_t *cblocks;
	int n_comp_blocks;
	int cblock_idx;
	double maf_limit;
	pthread_mutex_t cblock_mut;
	pthread_cond_t cblock_cond[3];
	pthread_mutex_t param_mut;
	pthread_cond_t param_cond;
	pthread_mutex_t contig_queue_mut;
	pthread_cond_t contig_queue_cond;
	int n_unsorted_being_processed;
	uint16_t n_prefix;
	bool input_finished;
	bool output_finished;
	bool unsorted_flag;
	bool abort_flag;
} dbsnp_param_t;

void handle_command_line(int argc, char *argv[], dbsnp_param_t * const par);
bool open_outfile(dbsnp_param_t *par);
void finish_output(dbsnp_param_t * const par);
void *input_thread(void *pt);
void *output_thread(void *pt);
void *compress_thread(void *pt);
void *write_thread(void *pt);
int add_to_bin(bin * const b, const snp_t * const snp, const int pref_ix, dbsnp_param_t * const par);
void clear_bins(bin *b, int ct);
void check_prefix(prefix **ppt, char *pref, int len, dbsnp_param_t * const par);
void add_remaining_contigs_to_queue(dbsnp_param_t * const par);
void free_used_files(dbsnp_param_t * const par);
void init_buffer(buffer_t * const buf, size_t sz);
buffer_t *new_buffer(size_t sz);
void destroy_buffer(buffer_t * const buf);
void resize_buffer(buffer_t * const buf, size_t sz);
void swap_buffers(buffer_t * const buf1, buffer_t * const buf2);
void read_alias_file(dbsnp_param_t * const par);
char *check_alias(snp_t * const snp, dbsnp_param_t * const par);
void read_select_file(dbsnp_param_t * const par);

#define reserve_buffer(buf, x) resize_buffer(buf, (buf)->len + x)
#define clear_buffer(buf) (buf)->len = 0

#endif /* dbSNP_idx_h */
