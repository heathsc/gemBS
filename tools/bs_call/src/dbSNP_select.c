/*
 * dbSNP_select.c
 *
 *  Created on: Feb 14, 2020
 *      Author: heath
 */

#include "dbSNP_idx.h"
#include "dbSNP_utils.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <ctype.h>

#include <htslib/khash.h>

KHASH_SET_INIT_STR(str);

void read_select_file(dbsnp_param_t * const par) {
	bool filter;
	FILE *fp = open_readfile_and_check(par->select_file, &filter);
	char *buf = NULL;
	size_t buf_size = 0;
	int n_snps = 0;
	khash_t(str) *h;
	h = kh_init(str);
	fprintf(stderr, "Reading selected SNP list from %s\n", par->select_file);
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
	par->select_hash = h;
	fclose(fp);
	if(buf) free(buf);
	if(filter) {
		int i;
		while(waitpid(-1, &i, WNOHANG) > 0);
	}
	fprintf(stderr, "List of %d unique SNPs read in\n", n_snps);
}
