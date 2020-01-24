/*
 * rec.c
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
#include <math.h>
#include <stdbool.h>
#include <pthread.h>
#include <time.h>
#include <sys/wait.h>
#include "utils.h"
#include "mextr.h"

rec_t *rec_init(const int ns) {
	assert(ns > 0);
	rec_t *rec = calloc((size_t)1, sizeof(rec_t));
	rec->sample_gt = malloc(sizeof(gt_meth) * ns);
	return rec;
}

#define REC_BUF_SHIFT 64

void *handle_rec_buf(void *p) {
	args_t * const args = p;
	rec_buffer_t *rec_buf = &args->rec_buf;
	rec_t * lbuf[REC_BUF_SIZE];
	const struct timespec wt = {0, 2500};
	int ix = 0;
	while(!args->proc_finished) {
		while(ix < REC_BUF_SIZE && ((rec_buf->buf[ix]->tasks & (REC_READY | RJ_ALL)) == REC_READY)) ix++;
		if(ix >= REC_BUF_SHIFT) {
			pthread_mutex_lock(&rec_buf->mut);
			rec_t ** bp = rec_buf->buf;
			for(int i = 0; i < ix; i++) {
				lbuf[i] = bp[i];
				lbuf[i]->tasks = 0;
			}
			for(int i = 0; i < REC_BUF_SIZE - ix; i++, bp++) *bp = bp[ix];
			for(int i = 0; i < ix; i++) *bp++ = lbuf[i];
			rec_buf->first_index += ix;
			pthread_mutex_unlock(&rec_buf->mut);
			ix = 0;
		} else nanosleep(&wt, NULL);
	}
	return NULL;
}

