/*
 * init_parampar->c
 *
 *  Created on: Nov 24, 2019
 *      Author: heath
 */

#include <stdio.h>

#include "gem_tools.h"
#include "bs_call.h"

void init_param(sr_param * const par) {

	const char *flts[4]= { "q20", "qd2", "fs60", "mq40" };
	const bool het[10] = { false, true, true, true, false, true, true, false, true, false };

	memset(par, 0, sizeof(sr_param));
	par->mmap_input = false;
	par->verbose = false;
	par->haploid = false;
	par->blank_trim = false;
	par->keep_unmatched = false;
	par->caller = maximum_likelihood;
	par->out_file_type = FT_UNKN;
	par->mapq_thresh = DEFAULT_MAPQ_THRESH;
	par->min_qual = MIN_QUAL;
	par->max_template_len = DEFAULT_MAX_TEMPLATE_LEN;
	par->under_conv = DEFAULT_UNDER_CONVERSION;
	par->over_conv = DEFAULT_OVER_CONVERSION;
	par->ref_bias = DEFAULT_REF_BIAS;
	par->keep_duplicates = false;
	par->ignore_duplicates = false;
	par->all_positions = false;
	par->benchmark_mode = false;
	par->work.print_end = par->work.process_end = false;
	par->work.vcf_ctg = NULL;
	pthread_mutex_init(&par->work.print_mutex, NULL);
	pthread_cond_init(&par->work.print_cond1, NULL);
	pthread_cond_init(&par->work.print_cond2, NULL);
	pthread_mutex_init(&par->work.vcf_mutex, NULL);
	pthread_cond_init(&par->work.vcf_cond, NULL);
	pthread_mutex_init(&par->work.process_mutex, NULL);
	pthread_cond_init(&par->work.process_cond1, NULL);
	pthread_cond_init(&par->work.process_cond2, NULL);
	pthread_mutex_init(&par->work.mprof_mutex, NULL);
	pthread_cond_init(&par->work.mprof_cond1, NULL);
	pthread_cond_init(&par->work.mprof_cond2, NULL);
	pthread_mutex_init(&par->work.calc_mutex, NULL);
	pthread_cond_init(&par->work.calc_cond1, NULL);
	pthread_cond_init(&par->work.calc_cond2, NULL);
	defs_t * const defs = &par->defs;
	for(int i = 0; i < 4 ;i++) defs->flt_name[i] = strdup(flts[i]);
	for(int i = 0; i < 10; i++) defs->gt_het[i] = het[i];
	lfact_store_init(defs->lfact_store);
	for(int i = 0; i < 100; i++) defs->logp[i] = log(0.01 * (double)(i + 1));
	uint8_t *tab = par->work.flt_tab;
	for(int q = MIN_QUAL; q < FLT_QUAL; q++) {
		int x = q << 2;
		// Strand 0 (non converted)
		tab[x] = 11; tab[x + 1] = 6; tab[x + 2] = 10; tab[x + 3] = 7;
		// Strand 1 (C2T)
		x += 256;
		tab[x] = 11; tab[x + 1] = 4; tab[x + 2] = 10; tab[x + 3] = 5;
		// Strand 2 (G2A)
		x += 256;
		tab[x] = 9; tab[x + 1] = 6; tab[x + 2] = 8; tab[x + 3] = 7;
	}
};
