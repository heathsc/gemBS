/*
 * unpack.c
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

#include "htslib/hfile.h"
#include "htslib/bgzf.h"

static void process(bcf1_t * const rec, rec_t * const rc, char ** cxp, int32_t * cx_np, const bool silent, args_t * const args)
{
	static char *ftags[N_FTAGS] = {"FT", "MC8", "AMQ", "CX", "AQ", "MQ"};
	static int ftypes[N_FTAGS] = {BCF_HT_STR, BCF_HT_INT, BCF_HT_INT, BCF_HT_STR, BCF_HT_INT, BCF_HT_INT};

	static int old_id = -1;
	static int old_complete = 0;

	int ns = bcf_hdr_nsamples(args->hdr);
	bcf_unpack(rec, BCF_UN_FLT);
	int n_all = rec->n_allele;
	bool cg = false;
	for(int i = 0; i < n_all; i++) {
		char c = rec->d.allele[i][0];
		if((c == 'C' || c == 'G') && rec->d.allele[i][1] == 0) {
			cg = true;
			break;
		}
	}
	rc->tasks = args->job_mask | REC_SKIP;
	if(cg) { // Site with potentially Cs or Gs
		bcf_unpack(rec, BCF_UN_ALL);
		// Get format tags
		for(int ix = 0; ix < N_FTAGS; ix++) {
			fmt_store_t *s = rc->tags + ix;
			s->ne = bcf_get_format_values(args->hdr, rec, ftags[ix], &s->dat_p, &s->dat_n, ftypes[ix]);
		}
		if(rc->tags[FMT_CX].ne > 0 && rc->tags[FMT_MC8].ne == ns * 8) {
			// Get sample base counts and genotype probs.
			int32_t *mc8_p = rc->tags[FMT_MC8].dat_p;
			int32_t *amq_p = rc->tags[FMT_AMQ].dat_p;
			int n_amq = rc->tags[FMT_AMQ].ne / ns;
			int32_t *aq_p = rc->tags[FMT_AQ].ne == ns ? rc->tags[FMT_AQ].dat_p : NULL;
			int32_t *mq_p = rc->tags[FMT_MQ].ne == ns ? rc->tags[FMT_MQ].dat_p : NULL;
			double ms_mq = 0.0;
			int32_t tot_n = 0;
			for(int i = 0; i < ns; i++) {
				int32_t *ct = rc->sample_gt[i].counts;
				int32_t *amq = rc->sample_gt[i].aqual;
				memset(ct, 0, sizeof(int32_t) * 8);
				memset(amq, 0, sizeof(int32_t) * 8);
				int32_t x = mc8_p[i * 8];
				int k = 0;
				if(x != bcf_int32_missing) {
					int k1 = 0;
					for(int j = 0; j < 8; j++) {
						x = mc8_p[i * 8 + j];
						ct[j] += x;
						k += x;
						if(x > 0 && amq_p != NULL && k1 < n_amq) {
							int q = amq_p[i * n_amq + k1++];
							if(q >= 0) {
								if(q > MAX_QUAL) q = MAX_QUAL;
								amq[j] = q;
							}
						}
					}
					if(amq_p == NULL) {
						int q = aq_p == NULL ? args->bq_thresh : aq_p[i];
						if(q > MAX_QUAL) q = MAX_QUAL;
						for(int j = 0; j < 8; j++) amq[j] = q;
					}
				}
				if(k > 0) {
					if(mq_p != NULL) {
						int m = mq_p[i];
						ms_mq += (double)k * (double)(m * m);
					}
					tot_n += k;
					calc_gt_prob(rc->sample_gt + i, args, rec->d.allele[0][0]);
					rc->sample_gt[i].skip = false;
				} else rc->sample_gt[i].skip = true;
			}
			// If we force a common genotype, calculate prob. distribution for common genotype
			if(ns > 1) {
				double gt[10];
				for(int k = 0; k < 10; k++) gt[k] = 0.0;
				for(int i = 0; i < ns; i++) {
					if(!rc->sample_gt[i].skip) {
						for(int k = 0; k < 10; k++) gt[k] += rc->sample_gt[i].gt_prob[k];
					}
				}
				double max = gt[0];
				int max_gt = 0;
				for(int k = 1; k < 10; k++) {
					if(gt[k] > max) {
						max = gt[k];
						max_gt = k;
					}
				}
				double sum = 0.0;
				for(int k = 0; k < 10; k++) sum += exp(gt[k] - max);
				sum = log(sum);
				for(int k = 0; k < 10; k++) gt[k] -= (max + sum);
				if(args->common_gt) {
					for(int i = 0; i < ns; i++) {
						if(!rc->sample_gt[i].skip) {
							for(int k = 0; k < 10; k++) rc->sample_gt[i].gt_prob[k] = gt[k];
							rc->sample_gt[i].max_gt = max_gt;
							rc->sample_gt[i].sum = max + sum;
						}
					}
				}
				rc->max_common_gt = max_gt;
			} else rc->max_common_gt = rc->sample_gt->max_gt;

			int cx_len = bcf_get_info_values(args->hdr, rec, "CX", (void **)cxp, cx_np, BCF_HT_STR);
			rc->ref = (cx_len >= 3 ? (*cxp)[2] : '.');
			int i;
			for(i = 0; i < cx_len && i < 5; i++) rc->cx[i] = (*cxp)[i];
			for(;i < 5; i++) rc->cx[i] = 'N';
			rc->pos = rec->pos;
			rc->rid = rec->rid;
			int id = args->id_trans[rc->rid];
			if(!silent) {
				const bool tty = isatty(STDERR_FILENO);
				int complete = 0;
				if(tty) {
					const uint64_t tot1 = rec->pos + 1 + (id ? args->cumul_len[id - 1] : 0);
					const uint64_t tot2 = args->cumul_len[args->sr->regions->nseqs - 1];
					complete = (int)((tot1 * 1000) / tot2);
				}
				if(id > old_id || complete > old_complete) {
					old_id = id;
					old_complete = complete;
					assert(id >= 0);
					if(tty) {
						fprintf(stderr,"Reading %s (%.1f%% completed)                                                      \r", args->sr->regions->seq_names[id], 0.1 * (double)complete);
					} else {
						fprintf(stderr,"Reading %s\n", args->sr->regions->seq_names[id], complete);
					}
				}
			}
			rc->tasks &= ~REC_SKIP;
		}
	}
	rc->tasks |= REC_READY;
}

#define BUF_SPACE(r, w) ( (w) >= (r) ? (r) + READ_BUF_SIZE - (w) : (r) - (w) )
#define NUM_READY(r, w) ( (w) >= (r) ? (w) - (r) : (w) + READ_BUF_SIZE - (r) )

#define READ_BLK_SIZE 64
#define UNPACK_BLK_SIZE 32

void *unpack_bcf_thread(void *p)
{
	gthr_info_t * const gi = p;
	args_t * const args = gi->args;
	const int thr_idx = gi->thread_idx;
	bcf1_t *lbuf[UNPACK_BLK_SIZE];
	rec_t * lrbuf[UNPACK_BLK_SIZE];
	uint64_t idx[UNPACK_BLK_SIZE];
	for(int i = 0; i < UNPACK_BLK_SIZE; i++) lbuf[i] = bcf_init();
	bcf1_buffer_t * const rb = &args->read_buf;
	rec_buffer_t * const rec_buf = &args->rec_buf;
	volatile int * const wp = &rb->write_pos;
	volatile int * const rp = &rb->read_pos;
	const struct timespec wt = {0, 5000};
	char *cx = NULL;
	int32_t cx_n = 0;

	while(true) {
		pthread_mutex_lock(&rb->mut);
		int k;
		while((k = NUM_READY(*rp, *wp)) < UNPACK_BLK_SIZE && !args->input_finished) {
			pthread_cond_wait(&rb->cond[1], &rb->mut);
		}
		if(!k) {
			pthread_mutex_unlock(&rb->mut);
			break;
		}
		if(k > UNPACK_BLK_SIZE) k = UNPACK_BLK_SIZE;
		for(int i = 0; i < k; i++) {
			bcf1_t *tmp = lbuf[i];
			lbuf[i] = rb->buf[*rp];
			idx[i] = rb->idx[*rp];
			rb->buf[*rp] = tmp;
			*rp = ((*rp) + 1) % READ_BUF_SIZE;
		}
		pthread_mutex_unlock(&rb->mut);
		pthread_cond_signal(&rb->cond[0]);
		// Check if we have space to store in rec_buf
		const uint64_t last_idx = idx[k - 1];
		pthread_mutex_lock(&rec_buf->mut);
		while(last_idx - rec_buf->first_index >= REC_BUF_SIZE ) {
			pthread_mutex_unlock(&rec_buf->mut);
			nanosleep(&wt, NULL);
			pthread_mutex_lock(&rec_buf->mut);
		}
		int ix = (int)(idx[0] - rec_buf->first_index);
		for(int i = 0; i < k; i++) lrbuf[i] = rec_buf->buf[ix + i];
		pthread_mutex_unlock(&rec_buf->mut);
		for(int i = 0; i < k; i++) {
			process(lbuf[i], lrbuf[i], &cx, &cx_n, thr_idx > 0, args);
		}
	}
	for(int i = 0; i < UNPACK_BLK_SIZE; i++) bcf_destroy(lbuf[i]);
	if(!thr_idx) {
		if(isatty(STDERR_FILENO)) fprintf(stderr,"Reading 100%% completed                                               \n");
		else fprintf(stderr, "Input finished\n");
	}
	return NULL;
}

void *read_thread(void * const p) {
	args_t * const args = p;
	bcf_srs_t * const sr = args->sr;
	bcf1_buffer_t * const rb = &args->read_buf;
	volatile int * const rp = &rb->read_pos;
	const struct timespec wt = {0, 5000};
	uint64_t count = 0;
	for(int i = 0; i < READ_BUF_SIZE; i++) rb->buf[i] = bcf_init();
	while(!args->input_finished) {
		int ix = rb->write_pos;
		pthread_mutex_lock(&rb->mut);
		while(BUF_SPACE(*rp, ix) < READ_BLK_SIZE) {
			pthread_cond_wait(&rb->cond[0], &rb->mut);
//			nanosleep(&wt, NULL);
		}
		int ix_end = (*rp);
		pthread_mutex_unlock(&rb->mut);
		int ix1 = (ix + 1) % READ_BUF_SIZE;
		while(ix1 != ix_end) {
			if(bcf_sr_next_line(sr)) {
				rb->idx[ix] = count++;
				bcf_sr_swap_line(sr, 0, rb->buf[ix]);
			} else {
				args->input_finished = true;
				break;
			}
			rb->write_pos = ix = ix1;
			ix1 = (ix + 1) % READ_BUF_SIZE;
		}
		pthread_cond_broadcast(&rb->cond[1]);
	}
	int ix = rb->write_pos;
	pthread_mutex_lock(&rb->mut);
	while((*rp) != ix) {
		pthread_cond_broadcast(&rb->cond[1]);
		pthread_cond_wait(&rb->cond[0], &rb->mut);
//		pthread_mutex_unlock(&rb->mut);
//		nanosleep(&wt, NULL);
//		pthread_mutex_lock(&rb->mut);
	}
	pthread_mutex_unlock(&rb->mut);
	for(int i = 0; i < READ_BUF_SIZE; i++) bcf_destroy(rb->buf[i]);
	return NULL;
}


