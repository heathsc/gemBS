/*
 * init_params.c
 *
 *  Created on: Dec 26, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <stdbool.h>
#include <pthread.h>
#include "utils.h"
#include "mextr.h"
#include "bbi.h"

#include "htslib/hfile.h"

void init_params(args_t *const args) {
	memset(args, 0, sizeof(args_t));
	args->bedmethyl_desc = ".";
	args->stats = NULL;
	args->min_prop = 0.0;
	args->min_num = 1;
	args->min_inform = 0;
	args->min_nc = 1;
	args->threads = 1;
	args->ref_bias = DEFAULT_REF_BIAS;
	args->under_conv = DEFAULT_UNDER_CONV;
	args->over_conv = DEFAULT_OVER_CONV;
	args->bq_thresh = DEFAULT_BQ_THRESH;
	args->mq_thresh = DEFAULT_MAPQ_THRESH;
	args->mode = CPGMODE_COMBINED;
	args->sel_mode = SELECT_HOM;
	args->sel_thresh = DEFAULT_SELECT_THRESH;
	args->header = true;
	pthread_mutex_init(&args->read_buf.mut, NULL);
	pthread_mutex_init(&args->rec_buf.mut, NULL);
	pthread_mutex_init(&args->cblock_buf.mut, NULL);
	for(int i = 0; i < 3; i++) {
		pthread_cond_init(args->cblock_buf.cond + i, NULL);
	}
	for(int i = 0; i < 2; i++) pthread_cond_init(args->read_buf.cond + i, NULL);
	for(int i = 0; i < 5; i++) {
		args->bb_global[i].buffer = malloc(sizeof(kstring_t));
		ks_initialize(args->bb_global[i].buffer);
		args->bb_global[i].first_time = true;
		args->bb_global[i].zoom_scales[0] = i < 3 ? INITIAL_REDUCTION : BW_INITIAL_REDUCTION;
		for(int j = 1; j < ZOOM_LEVELS; j++) args->bb_global[i].zoom_scales[j] = args->bb_global[i].zoom_scales[j - 1] * ZOOM_RES_INCREMENT;
	}
	for(int i = 0; i < 3; i++) ks_initialize(args->pr_str + i);
}
