/*
 * dbSNP_input.c
 *
 *  Created on: Feb 2, 2020
 *      Author: heath
 */

#include "dbSNP_idx.h"
#include "dbSNP_utils.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <ctype.h>
#include "../include/dbSNP_json.h"

static file_t *next_file(dbsnp_param_t * const par) {
	if(par->abort_flag) return NULL;
	pthread_mutex_lock(&par->param_mut);
	file_t *file = par->unsorted;
	if(file) HASH_DEL(par->unsorted, file);
	else {
		while(par->unsorted_flag && par->sorted) pthread_cond_wait(&par->param_cond, &par->param_mut);
		file = par->sorted;
		if(file) HASH_DEL(par->sorted, file);
	}
	if(file) {
		if(!file->sorted) par->n_unsorted_being_processed++;
		file->hh.next = par->used_files;
		par->used_files = file;
		pthread_mutex_unlock(&par->param_mut);
		if(strcmp(file->name, "-")) file->fp = open_readfile(file->name, &file->filter);
		else {
			file->fp = stdin;
			file->filter = false;
		}
	} else 	pthread_mutex_unlock(&par->param_mut);
	return file;
}

void free_used_files(dbsnp_param_t * const par) {
	file_t *f = par->used_files, *tmp;
	while(f) {
		tmp = f->hh.next;
		free(f);
		f = tmp;
	}
}

static void set_header(dbsnp_param_t * const par, const char * const hd) {
	pthread_mutex_lock(&par->param_mut);
	if(!par->header) par->header = strdup(hd);
	pthread_mutex_unlock(&par->param_mut);
}

static contig *new_contig(const snp_t * const snp, char * const rname,  file_t * const file, const uint32_t binx) {
	contig * const ctg = malloc(sizeof(contig));
	ctg->min_bin = ctg->max_bin = binx;
	ctg->bins = malloc(sizeof(bin));
	ctg->name = malloc(snp->name_len + 1);
	ctg->name_len = snp->cname_len;
	memcpy(ctg->name, snp->cname, ctg->name_len);
	ctg->name[ctg->name_len] = 0;
	ctg->rname = rname ? rname : ctg->name;
	ctg->first_file = file;
	ctg->in_queue = false;
	pthread_mutex_init(&ctg->mut, NULL);
	clear_bins(ctg->bins, 1);
	return ctg;
}

static contig * find_or_create_contig(dbsnp_param_t * const par, const snp_t * const snp, char * const rname, file_t * const file, const uint32_t binx) {
	contig *ctg;
	pthread_mutex_lock(&par->param_mut);
	HASH_FIND(hh, par->contigs, snp->cname, snp->cname_len, ctg);
	if(!ctg) {
		ctg = new_contig(snp, rname, file, binx);
		HASH_ADD_KEYPTR(hh, par->contigs, ctg->name, ctg->name_len, ctg);
	} else if(!par->unsorted_flag) {
		if(file == ctg->first_file) {
			fprintf(stderr, "Contig %s found in multiple places in input '%s' - verify usage of --sorted option\n", ctg->rname, file->name);
			ctg = NULL;
		} else {
			file_t *file1= ctg->first_file;
			if(file1->sorted || !file1->read) {
				fprintf(stderr, "Contig %s found in multiple files ('%s' and '%s') - verify usage of --sorted option\n", ctg->rname, file->name, ctg->first_file->name);
				ctg = NULL;
			} else ctg->first_file = file;
		}
	}
	pthread_mutex_unlock(&par->param_mut);
	return ctg;
}

static void add_contig_to_queue(dbsnp_param_t * const par, contig * const ctg) {
	if(!par->unsorted_flag) {
		pthread_mutex_lock(&par->contig_queue_mut);
		ctg->next = par->contig_queue;
		par->contig_queue = ctg;
		ctg->in_queue = true;
		pthread_cond_signal(&par->contig_queue_cond);
		pthread_mutex_unlock(&par->contig_queue_mut);
	}
}

void add_remaining_contigs_to_queue(dbsnp_param_t * const par) {
	for(contig *ctg = par->contigs; ctg; ctg = ctg->hh.next) {
		if(!ctg->in_queue) add_contig_to_queue(par, ctg);
	}
}

dbsnp_input_type_t guess_input_type(char * const buf, const ssize_t l) {
	dbsnp_input_type_t itype = dbsnp_bed;
	if(buf[0] == '{') itype = dbsnp_json;
	else if(l >= 16 && !strncmp(buf, "##fileformat=VCF",16)) itype = dbsnp_vcf;
	return itype;
}

static tokens *parse_bed_line(char * const buf, const ssize_t l, tokens *tok, snp_t * const snp, dbsnp_param_t * const par) {
	if(l > 6 && !strncmp(buf, "track ", 6)) {
		if(par->header == NULL) set_header(par, buf);
	} else {
		tok = tokenize(buf, '\t', tok);
		if(tok->n_tok > 4) {
			char *p;
			uint32_t x = (uint32_t)strtoul(tok->toks[1], &p, 10);
			uint32_t y = (uint32_t)strtoul(tok->toks[2], &p, 10);
			if(y > x && y - x == 1) {
				snp->pos = y;
				snp->cname = tok->toks[0];
				snp->cname_len = tok->toks[1] - tok->toks[0] - 1;
				snp->name = tok->toks[3];
				snp->name_len = tok->toks[4] - tok->toks[3] - 1;
				snp->ok = true;
			}
		}
	}
	return tok;
}

static tokens *parse_vcf_line(char * const buf, const ssize_t l, tokens *tok, snp_t * const snp, dbsnp_param_t * const par) {
	if(buf[0] != '#') {
		tok = tokenize(buf, '\t', tok);
		if(tok->n_tok > 4 && tok->toks[3][1] == 0 && tok->toks[4][1] == 0) {
			char *p;
			snp->pos = (uint32_t)strtoul(tok->toks[1], &p, 10);
			snp->cname = tok->toks[0];
			snp->cname_len = tok->toks[1] - tok->toks[0] - 1;
			snp->name = tok->toks[2];
			snp->name_len = tok->toks[3] - tok->toks[2] - 1;
			snp->ok = true;
		}
	}
	return tok;
}

static void adjust_name(snp_t * const snp, prefix **pref, dbsnp_param_t * const par) {
	int k = snp->name_len;
	for(; k > 0; k--) if(!isdigit((int)snp->name[k - 1])) break;
	check_prefix(pref, snp->name, k, par);
	snp->name_off = k;
}

void *input_thread(void *pt) {
	dbsnp_param_t * const par = pt;
	bool *st = malloc(sizeof(bool));
	*st = false;
	char *buf = NULL;
	size_t buf_size = 0;
	tokens *tok = NULL;
	contig *ctg = NULL;
	prefix *pref = NULL;
	jsmn_work_t jwork = {
			.keys = NULL,
			.tcount = 0,
			.tok = NULL
	};
	uint64_t n_snps = 0;
	while(!*st) {
		file_t * const file = next_file(par);
		if(!file) break;
		FILE *fin = file->fp;
		fprintf(stderr, "Reading from %s\n", fin != stdin ? file->name : "<stdin>");
		dbsnp_input_type_t itype = par->input_type;
		snp_t snp;
		while(!*st) {
			if(par->abort_flag) {
				*st = true;
				break;
			}
			ssize_t l = getline(&buf, &buf_size, fin);
			if(l < 0) break;
			if(l > 0 && buf[l - 1] == '\n') buf[--l] = 0;
			if(l > 0) {
				if(itype == dbsnp_auto) {
					itype = guess_input_type(buf, l);
					if(itype == dbsnp_json) check_prefix(&pref, "rs", 2, par);
				}
				snp.ok = false;
				snp.maf = -1.0;
				snp.name_off = 0;
				switch(itype) {
				case dbsnp_json:
					parse_json_line(buf, l, &jwork, &snp, par);
					break;
				case dbsnp_bed:
					tok = parse_bed_line(buf, l, tok, &snp, par);
					if(snp.ok) adjust_name(&snp, &pref, par);
					break;
				case dbsnp_vcf:
					tok = parse_vcf_line(buf, l, tok, &snp, par);
					if(snp.ok) adjust_name(&snp, &pref, par);
					break;
				default:
					fprintf(stderr, "Input type nor currently handled\n");
					*st = true;
					par->abort_flag = true;
					break;
				}
				if(*st) break;
				if(snp.ok) {
					uint32_t binx = snp.pos >> 6;
					if(ctg == NULL || ctg->name_len != snp.cname_len || memcmp(ctg->name, snp.cname, snp.cname_len)) {
						char *rname = NULL;
						if(par->aliases) {
							rname = check_alias(&snp, par);
							if(!rname) snp.ok = false;
						}
						if(snp.ok) {
							if(ctg) add_contig_to_queue(par, ctg);
							ctg = find_or_create_contig(par, &snp, rname, file, binx);
							if(!ctg) {
								*st = true;
								par->abort_flag = true;
								break;
							}
						}
					}
					if(snp.ok) {
						pthread_mutex_lock(&ctg->mut);
						if(binx > ctg->max_bin) {
							ctg->bins = realloc(ctg->bins, sizeof(bin) * (size_t)(binx - ctg->min_bin + 1));
							clear_bins(ctg->bins + (ctg->max_bin - ctg->min_bin + 1), binx - ctg->max_bin);
							ctg->max_bin = binx;
						} else if(binx < ctg->min_bin) {
							ctg->bins = realloc(ctg->bins, sizeof(bin) * (size_t)(ctg->max_bin - binx + 1));
							memmove(ctg->bins + (ctg->min_bin - binx), ctg->bins, sizeof(bin) * (size_t)(ctg->max_bin - ctg->min_bin + 1));
							clear_bins(ctg->bins, ctg->min_bin - binx);
							ctg->min_bin = binx;
						}
						n_snps += add_to_bin(ctg->bins + (binx - ctg->min_bin), &snp, pref->ix, par);
						pthread_mutex_unlock(&ctg->mut);
					}
				}
			}
		}
		if(fin != stdin) {
			fclose(fin);
			if(file->filter) {
				int i;
				while (waitpid(-1, &i, WNOHANG) > 0);
			}
		}
		file->read = true;
		if(!file->sorted) {
			pthread_mutex_lock(&par->param_mut);
			assert(par->unsorted_flag &&& par->n_unsorted_being_processed > 0);
			if(!(--par->n_unsorted_being_processed) && !par->unsorted) {
				par->unsorted_flag = false;
				fprintf(stderr, "Switching to sorted mode (ctg = %s)\n", ctg->name);
				pthread_cond_broadcast(&par->param_cond);
			}
			pthread_mutex_unlock(&par->param_mut);
			ctg = NULL;
		}
	}
	if(ctg) add_contig_to_queue(par, ctg);
	if(buf != NULL) free(buf);
	if(tok != NULL) free_tokens(tok);
	pthread_mutex_lock(&par->param_mut);
	par->n_snps += n_snps;
	pthread_mutex_unlock(&par->param_mut);
	return st;
}
