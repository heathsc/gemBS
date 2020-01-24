/*
 * snplist.c
 *
 *  Created on: Jan 23, 2020
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
#include <ctype.h>
#include <sys/wait.h>

#include "snpxtr.h"

#include <htslib/khash.h>

#include "utils.h"

KHASH_SET_INIT_STR(str);

void read_snplist(sargs_t * const args) {
	bool filter;
	FILE *fp = open_readfile_and_check(args->snplistname, &filter);
	char *buf = NULL;
	size_t buf_size = 0;
	int n_snps = 0;
	khash_t(str) *h;
	h = kh_init(str);
	fprintf(stderr, "Reading SNP list from %s\n", args->snplistname);
	for(;;) {
		ssize_t l = getline(&buf, &buf_size, fp);
		if(l < 0) break;
		int i;
		for(i = 0; i < l; i++) if(!isspace(buf[i])) break;
		int j = i;
		for(; i < l; i++) if(isspace(buf[i])) break;
		if(i == j) continue;
		buf[i] = 0;
		int not_found;
		khint_t k = kh_put(str, h, buf + j, &not_found);
		if(not_found) {
			n_snps++;
			kh_key(h, k) = strdup(buf + j);
		}
	}
	args->snp_hash = h;
	fclose(fp);
	if(buf) free(buf);
	if(filter) {
		int i;
		while(waitpid(-1, &i, WNOHANG) > 0);
	}
	fprintf(stderr, "List of %d unique SNPs read in\n", n_snps);
}

