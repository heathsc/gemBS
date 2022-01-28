//
//  dbSNP_idx.c
//  dbSNP
//
//  Created by Simon Heath on 15/11/2017.
//  Copyright 2017 Simon Heath. All rights reserved.
//

#include "dbSNP_idx.h"
#include "dbSNP_utils.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <zlib.h>

int main(int argc, char *argv[]) {
	dbsnp_param_t par = {
			.output_file = NULL,
			.outfile = NULL,
			.threads = 0,
			.input_type = dbsnp_auto,
			.header = NULL,
			.contigs = NULL,
			.contig_queue = NULL,
			.aliases = NULL,
			.prefixes = NULL,
			.n_prefix = 0,
			.n_snps = 0,
			.n_snps_maf_filtered = 0,
			.max_buf_size = 0,
			.maf_limit = 1.0,
			.param_mut = PTHREAD_MUTEX_INITIALIZER,
			.param_cond = PTHREAD_COND_INITIALIZER,
			.contig_queue_mut = PTHREAD_MUTEX_INITIALIZER,
			.contig_queue_cond = PTHREAD_COND_INITIALIZER,
			.cblock_mut = PTHREAD_MUTEX_INITIALIZER,
			.cblock_cond = {PTHREAD_COND_INITIALIZER, PTHREAD_COND_INITIALIZER, PTHREAD_COND_INITIALIZER},
			.cblock_idx = 0,
			.input_finished = false,
			.output_finished = false,
			.sorted = NULL,
			.unsorted = NULL,
			.select_file = NULL,
			.select_hash = NULL,
			.used_files = 0,
			.unsorted_flag = false,
			.n_unsorted_being_processed = 0
	};
	par.output_buf = new_buffer(4096);
	handle_command_line(argc, argv, &par);
	if(par.chrom_alias_file) read_alias_file(&par);
	bool st = open_outfile(&par);
	if(!st) {
		const int nt = par.threads = par.threads > 0 ? par.threads : 1;
		const int nf = HASH_COUNT(par.sorted) + HASH_COUNT(par.unsorted);
		par.unsorted_flag = HASH_COUNT(par.unsorted) > 0;
		const int read_jobs = nt > nf ? nf : nt;
		par.read_jobs = read_jobs;
		par.compress_jobs = nt;
		par.n_comp_blocks = par.compress_jobs << 3;
		par.cblocks = malloc(sizeof(comp_block_t) * par.n_comp_blocks);
		for(int i = 0; i < par.n_comp_blocks; i++) {
			par.cblocks[i].state = block_empty;
			init_buffer(&par.cblocks[i].cblock, 4096);
			init_buffer(&par.cblocks[i].ublock, 4096);
		}
		pthread_t write_thr;
		pthread_create(&write_thr, NULL, write_thread, &par);
		pthread_t *compress_thr = malloc(sizeof(pthread_t) * par.compress_jobs);
		for(int i = 0; i < par.compress_jobs; i++) pthread_create(compress_thr + i, NULL, compress_thread, &par);
		pthread_t output_thr;
		pthread_create(&output_thr, NULL, output_thread, &par);
		pthread_t *input_thr = malloc(sizeof(pthread_t) * read_jobs);
		for(int i = 0; i < read_jobs; i++) pthread_create(input_thr + i, NULL, input_thread, &par);
		bool *retval;
		for(int i = 0; i < read_jobs; i++) {
			pthread_join(input_thr[i], (void **)&retval);
			if(retval && *retval) st = true;
		}
		free(input_thr);
		assert(!par.unsorted_flag);
		free_used_files(&par);
		if(!st) add_remaining_contigs_to_queue(&par);
		pthread_mutex_lock(&par.contig_queue_mut);
		par.input_finished = true;
		pthread_cond_signal(&par.contig_queue_cond);
		pthread_mutex_unlock(&par.contig_queue_mut);
		pthread_join(output_thr, NULL);
		par.output_finished = true;
		pthread_cond_broadcast(&par.cblock_cond[1]);
		pthread_cond_broadcast(&par.cblock_cond[2]);
		for(int i = 0; i < par.compress_jobs; i++) pthread_join(compress_thr[i], NULL);
		free(compress_thr);
		pthread_join(write_thr, NULL);
		if(!st) finish_output(&par);
	}
	if(st) fprintf(stderr, "Index file creation failed\n");
	else {
		if(par.n_snps_maf_filtered > 0) fprintf(stderr, "Index file created: %lu snps processed, %lu selected for constant calling\n", par.n_snps, par.n_snps_maf_filtered);
		else fprintf(stderr, "Index file created: %lu snps processed\n", par.n_snps);
	}
	return st;
}
