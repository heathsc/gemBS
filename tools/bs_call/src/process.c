/*
 * process.c
 *
 *  Created on: Nov 24, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <sys/wait.h>
#include <pthread.h>
#include <time.h>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include "gem_tools.h"
#include "bs_call.h"

void *mprof_thread(void *arg) {
	sr_param * const par = arg;
	pthread_mutex_lock(&par->work.mprof_mutex);
	while(1) {
		while(par->work.mprof_read_idx == par->work.mprof_write_idx && !par->work.mprof_end) {
			struct timespec ts;
			clock_gettime(CLOCK_REALTIME, &ts);
			ts.tv_sec += 5;
			pthread_cond_timedwait(&par->work.mprof_cond1, &par->work.mprof_mutex, &ts);
		}
		bool end = (par->work.mprof_read_idx == par->work.mprof_write_idx);
		int ix = par->work.mprof_read_idx;
		pthread_mutex_unlock(&par->work.mprof_mutex);
		if(end) break;
		mprof_thread_t * mp = &par->work.mprof_thread[ix];
		meth_profile(mp->al, mp->x, mp->orig_pos, mp->max_pos, par);
		pthread_mutex_lock(&par->work.mprof_mutex);
		par->work.mprof_read_idx = (ix + 1) % N_MPROF_BUFFERS;
		pthread_cond_broadcast(&par->work.mprof_cond2);
	}
	return NULL;
}

void *process_thread(void *arg) {
	sr_param *par = arg;
	gt_vector *prev_align = NULL;
	work_t * const work = &par->work;
	while(true) {
		pthread_mutex_lock(&par->work.process_mutex);
		while(!par->work.align_list_waiting && !par->work.process_end) {
			struct timespec ts;
			clock_gettime(CLOCK_REALTIME, &ts);
			ts.tv_sec += 5;
			pthread_cond_timedwait(&par->work.process_cond1, &par->work.process_mutex, &ts);
		}
		pthread_mutex_unlock(&par->work.process_mutex);
		gt_vector *alist = work->align_list_waiting;
		if(alist != NULL) {
			ctg_t *ctg = work->ctg_waiting;
			uint64_t pos = work->y_waiting;
			work->ctg_waiting = NULL;
			work->free_list_waiting = prev_align;
			work->align_list_waiting = NULL;
			pthread_mutex_lock(&par->work.process_mutex);
			pthread_cond_signal(&work->process_cond2);
			pthread_mutex_unlock(&par->work.process_mutex);
			gt_status err = process_template_vector(alist, ctg, pos, par);
			if(err != GT_STATUS_OK) break;
			prev_align = alist;
		} else break;
	}
	return NULL;
}

void *print_thread(void *arg) {
	sr_param * const par = arg;
	bcf1_t *bcf = bcf_init();

	pthread_mutex_lock(&par->work.print_mutex);
	while(1) {
		while(!par->work.vcf_n && !par->work.print_end) {
			struct timespec ts;
			clock_gettime(CLOCK_REALTIME, &ts);
			ts.tv_sec += 5;
			pthread_cond_timedwait(&par->work.print_cond1, &par->work.print_mutex, &ts);
		}
		pthread_mutex_unlock(&par->work.print_mutex);
		if(par->work.vcf_n) {
			const char *ref_st = gt_string_get_string(par->work.ref);
			for(int i = 0; i < par->work.vcf_n; i++) {
				while(!par->work.vcf[i].ready) {
					pthread_mutex_lock(&par->work.vcf_mutex);
					while(!par->work.vcf[i].ready) {
						struct timespec ts;
						clock_gettime(CLOCK_REALTIME, &ts);
						ts.tv_sec += 5;
						pthread_cond_timedwait(&par->work.vcf_cond, &par->work.vcf_mutex, &ts);
					}
					pthread_mutex_unlock(&par->work.vcf_mutex);
				}
				print_vcf_entry(bcf, par->work.vcf_ctg, &par->work.vcf[i].gtm, ref_st, i + par->work.vcf_x,
						par->work.vcf_x, par->work.vcf[i].skip, par);
			}
			flush_vcf_entries(bcf, par);
			par->work.vcf_n = 0;
			pthread_mutex_lock(&par->work.print_mutex);
			pthread_cond_signal(&par->work.print_cond2);
		} else break;
	}
	return NULL;
}

gt_status open_input_file(sr_param * const param) {
	gt_status err = GT_STATUS_OK;
	htsFile *in_file = NULL;
	param->work.sam_idx = NULL;
	if(param->input_file != NULL) {
		in_file = hts_open(param->input_file, "r");
		if(in_file != NULL) {
			size_t l = strlen(param->name_reference_file);
			char *fai = gt_malloc(l + 5);
			sprintf(fai, "%s.fai", param->name_reference_file);
			int ret = hts_set_fai_filename(in_file, fai);
			if(ret) fprintf(stderr,"Could not open reference index %s\n", fai);
			free(fai);
			param->work.sam_idx = sam_index_load(in_file, param->input_file);
		}
	} else {
		hFILE *hf = hdopen(STDIN_FILENO, "r");
		if(hf != NULL) in_file = hts_hopen(hf, "<STDIN>", "r");
	}
	param->work.sam_file = in_file;
	if(!in_file) {
		err = GT_STATUS_FAIL;
		fprintf(stderr, "Could not open input file %s\n", param->input_file ? param->input_file : "<STDIN>");
	}
	return err;
}

void close_input_file(sr_param * const param) {
	if(param->work.sam_file) {
		hts_close(param->work.sam_file);
		if(param->work.sam_idx != NULL) hts_idx_destroy(param->work.sam_idx);
	}
}

gt_status bs_call_process(sr_param * const param) {
	gt_status err = GT_STATUS_OK;
	htsFile *in_file = param->work.sam_file;
	assert(in_file);
	bam_hdr_t *header = NULL;
	if(param->num_threads[INPUT_THREADS] > 0) hts_set_threads(in_file, param->num_threads[INPUT_THREADS]);
	fprintf(stderr,"Opened %s for input (%s)\n", in_file->fn, param->work.sam_idx ? "Index" : "No index");
	header = sam_hdr_read(in_file);
	if(header != NULL) {
		bool ret = process_sam_header(param, header);
		if(!ret) {
			print_vcf_header(param, header);
			param->work.ref = gt_string_new(16384);
			param->work.ref1 = gt_string_new(16384);
			for(int i = 0; i < N_MPROF_BUFFERS; i++) {
				param->work.mprof_thread[i].orig_pos[0] = gt_vector_new(256, sizeof(int));
				param->work.mprof_thread[i].orig_pos[1] = gt_vector_new(256, sizeof(int));
			}
			param->work.sam_header = header;
			fill_base_prob_table();
			init_calc_threads(param);
			pthread_t process_thr, print_thr, mprof_thr;
			pthread_create(&print_thr, NULL, print_thread, param);
			pthread_create(&process_thr, NULL, process_thread, param);
			pthread_create(&mprof_thr, NULL, mprof_thread, param);
			gt_vector *align_list = gt_vector_new(32, sizeof(align_details *));
			err = read_input(in_file, align_list, param);
			param->work.process_end = true;
			pthread_mutex_lock(&param->work.process_mutex);
			pthread_cond_broadcast(&param->work.process_cond1);
			pthread_mutex_unlock(&param->work.process_mutex);
			pthread_join(process_thr, 0);
			param->work.print_end = true;
			param->work.mprof_end = true;
			join_calc_threads(param);
			pthread_mutex_lock(&param->work.mprof_mutex);
			pthread_cond_broadcast(&param->work.mprof_cond1);
			pthread_cond_broadcast(&param->work.mprof_cond2);
			pthread_mutex_unlock(&param->work.mprof_mutex);
			pthread_mutex_lock(&param->work.print_mutex);
			pthread_cond_signal(&param->work.print_cond1);
			pthread_cond_signal(&param->work.print_cond2);
			pthread_mutex_unlock(&param->work.print_mutex);
			pthread_join(mprof_thr, 0);
			pthread_join(print_thr, 0);
			hts_close(param->work.vcf_file);
		}
	}
	close_input_file(param);
	return err;
}
