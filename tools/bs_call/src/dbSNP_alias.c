/*
 * dbSNP_alias.c
 *
 *  Created on: Feb 9, 2020
 *      Author: heath
 */

#include "dbSNP_idx.h"
#include "dbSNP_utils.h"
#include <sys/types.h>
#include <sys/wait.h>

void read_alias_file(dbsnp_param_t * const par) {
	bool filter;
	FILE *fp = open_readfile(par->chrom_alias_file, &filter);
	char *buf = NULL;
	size_t buf_size = 0;
	tokens *tok = NULL;
	while(true) {
		ssize_t l = getline(&buf, &buf_size, fp);
		if(l < 0) break;
		if(l > 0 && buf[l - 1] == '\n') buf[--l] = 0;
		if(l > 0) {
			tok = tokenize(buf, '\t', tok);
			if(tok->n_tok >= 2) {
				chrom_alias_t *alias;
				HASH_FIND_STR(par->aliases, tok->toks[0], alias);
				if(!alias) {
					alias = malloc(sizeof(chrom_alias_t));
					alias->alias = strdup(tok->toks[0]);
					alias->name = strdup(tok->toks[1]);
					HASH_ADD_KEYPTR(hh, par->aliases, alias->alias, strlen(alias->alias), alias);
				}

			}
		}
	}
	fclose(fp);
	if(filter) {
		int i;
		while (waitpid(-1, &i, WNOHANG) > 0);
	}
}

char *check_alias(snp_t * const snp, dbsnp_param_t * const par) {
	chrom_alias_t *alias = NULL;
	if(par->aliases) HASH_FIND(hh, par->aliases, snp->cname, snp->cname_len, alias);
	return alias ? alias->name : NULL;
}
