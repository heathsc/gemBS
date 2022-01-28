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
#include "bbi.h"

#include "htslib/hfile.h"
#include "htslib/bgzf.h"

void error(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
	exit(-1);
}

void destroy(args_t * const args)
{
	write_stats(args);
}

int main(int argc, char **argv) {

	args_t args;
	init_params(&args);
	handle_command_line(argc, argv, &args);
	int ns = bcf_hdr_nsamples(args.hdr);
	assert(ns > 0);
	init_files(&args);
	pthread_t cpg_thr, noncpg_thr, wig_thr, bedmethyl_thr;
	thr_info_t noncpg_info = {&args, output_noncpg, RJ_OUTPUT_NONCPG};
	thr_info_t bedmethyl_info = {&args, output_bedmethyl, RJ_OUTPUT_BEDMETHYL};
	if(args.cpgfile) {
		args.job_mask |= RJ_OUTPUT_CPG;
		pthread_create(&cpg_thr, NULL, cpg_thread, &args);
	}
	if(args.noncpgfile) {
		args.job_mask |= RJ_OUTPUT_NONCPG;
		pthread_create(&noncpg_thr, NULL, output_thread, &noncpg_info);
	}
	if(args.bedmethyl) {
		args.job_mask |= RJ_OUTPUT_BEDMETHYL;
		pthread_create(&bedmethyl_thr, NULL, handle_bedmethyl_thread, &bedmethyl_info);
	}

	if(args.header || args.bedmethyl) print_headers(&args);
	if(args.reportfilename != NULL) init_stats(&args);

	args.sample_Q[0] = malloc(sizeof(double) * ((ns + 1) * 2 + ns));
	args.sample_Q[1] = args.sample_Q[0] + ns + 1;
	args.sample_Q[2] = args.sample_Q[1] + ns + 1;
	args.sample_Q1[0] = malloc(sizeof(double) * ((ns + 1) * 2 + ns));
	args.sample_Q1[1] = args.sample_Q1[0] + ns + 1;
	args.sample_Q1[2] = args.sample_Q1[1] + ns + 1;
	args.sample_cpg = (args.mode == CPGMODE_COMBINED) ? malloc(sizeof(cpg_prob) * ns) : NULL;
	fill_base_prob_table();
	pthread_t read_thr;
	pthread_create(&read_thr, NULL, read_thread, &args);
	// If we are generating bigBed and bigWig files, we will take half of the threads for the compression stage
	if(args.bedmethyl && args.threads > 1) args.threads = (args.threads + 1) >> 1;
	int nt = 1;
	// If the user has asked for more threads we will take one extra thread for the bcf unpacking - more than this is rarely useful
	if(args.threads > 4) {
		nt++;
		args.threads--;
	}
	pthread_t * const unpack_bcf_thr = malloc(sizeof(pthread_t) * nt);
	gthr_info_t * const gthr_info = malloc(sizeof(gthr_info_t) * nt);
	for(int i = 0; i < nt; i++) {
		gthr_info[i].args = &args;
		gthr_info[i].thread_idx = i;
		pthread_create(unpack_bcf_thr + i, NULL, unpack_bcf_thread, gthr_info + i);
	}
	pthread_t handle_rec_buf_thr;
	pthread_create(&handle_rec_buf_thr, NULL, handle_rec_buf, &args);
	pthread_join(read_thr, NULL);
	for(int i = 0; i < nt; i++) pthread_join(unpack_bcf_thr[i], NULL);
	args.proc_finished = true;
	pthread_join(handle_rec_buf_thr, NULL);
	args.rec_finished = true;
	if(args.cpgfilename) pthread_join(cpg_thr, NULL);
	if(args.noncpgfilename) pthread_join(noncpg_thr, NULL);
	if(args.bedmethyl) pthread_join(bedmethyl_thr, NULL);
	fprintf(stderr,"mextr: finished\n");
	destroy(&args);
	bcf_sr_destroy(args.sr);
	return 0;
}


