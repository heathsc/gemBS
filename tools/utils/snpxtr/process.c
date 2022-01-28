/*
 * process.c
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
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <pthread.h>

#include "utils.h"
#include "snpxtr.h"

#include <htslib/hts.h>
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/tbx.h"
#include <htslib/khash.h>

static int base_tab[256] = {
		['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4
};
static void * dat_p = NULL;
static int dat_n = 0;

KHASH_SET_INIT_STR(str);

void make_tabix_index(char * const fname) {
	tbx_conf_t conf = tbx_conf_vcf;
	tbx_index_build(fname, 0, &conf);
}

void *md5_thread(void *p) {
	calc_file_md5(p);
	return NULL;
}

void process_input(sargs_t * const args) {
	dbsnp_header_t * const dbsnp_hdr = args->dbSNP_hdr;
	dbsnp_ctg_t *dbSNP_ctg = NULL;
	char rs[512];
	int rid = -1;
	bcf_srs_t * const sr = args->sr;
	bcf1_t *rec = bcf_init();
	khash_t(str) *h = args->snp_hash;

	while(bcf_sr_next_line(sr)) {
		bcf_sr_swap_line(sr, 0, rec);
		if(rec->rid != rid) {
			if(dbsnp_hdr) {
				const char * const cname = args->hdr->id[BCF_DT_CTG][rec->rid].key;
				if(dbSNP_ctg) unload_dbSNP_ctg(dbSNP_ctg);
				HASH_FIND(hh, dbsnp_hdr->dbSNP, cname, strlen(cname), dbSNP_ctg);
				if(dbSNP_ctg) {
					bool ret = load_dbSNP_ctg(dbsnp_hdr, dbSNP_ctg);
					if(!ret) dbSNP_ctg = NULL;
				}
			}
			rid = rec->rid;
		}
		int ns = bcf_hdr_nsamples(args->hdr);
		bcf_unpack(rec, BCF_UN_ALL);
		char *id = rec->d.id;
		bool rs_id;
		if(id[0]=='.' && id[1]==0) {
			rs_id = false;
			if(dbSNP_ctg) {
				if(dbSNP_lookup_name(dbsnp_hdr, dbSNP_ctg, rs, NULL, rec->pos + 1)) {
					id = rs;
					rs_id = true;
				}
			}
		} else rs_id = true;
		if(rs_id) {
			bool passed = true;
			if(h) {
				khint_t k = kh_get(str, h, id);
				if(k == kh_end(h)) passed = false;
			}
			if(passed) {
				passed = false;
				for(int i = 0; i < rec->d.n_flt; i++) if(rec->d.flt[i] == args->pass_idx) {
					passed = true;
					break;
				}
			}
			// Check for SNP
			int n_all = rec->n_allele;
			if(passed) {
				if(n_all > 4) passed = false;
				else for(int i = 0; i < n_all; i++) {
					char * const p = rec->d.allele[i];
					if(p[1] || !base_tab[(int)p[0]]) {
						passed = false;
						break;
					}
				}
			}
			if(passed) {
				int ne = bcf_get_format_values(args->hdr, rec, "FT", &dat_p, &dat_n, BCF_HT_STR);
				int gt_ix = -1;
				bcf_fmt_t * const fmt = rec->d.fmt;
				// Get GT Tag
				for(int i = 0; i < (int)rec->n_fmt; i++) {
					if(!fmt[i].p) continue;
					if(fmt[i].id == args->gt_idx) {
						gt_ix = i;
						break;
					}
				}
				if(gt_ix >= 0) {
					int sz = ne / ns;
					char *flt = dat_p;
					bcf_fmt_t *fmt = rec->d.fmt + gt_ix;
					passed = false;
					for(int i = 0; i < ns; i++) {
						args->gt[i] = 0;
						switch(fmt->type) {
						case BCF_BT_INT8:
						{
							if(fmt->n == 2) {
								int8_t *p = (int8_t *)(fmt->p + i * fmt->size);
								if(p[0] != bcf_int8_vector_end && p[1] != bcf_int8_vector_end) {
									int a1 = p[0] >> 1;
									int a2 = p[1] >> 1;
									if(a1 < 1 || a1 > n_all || a2 < 1 || a2 > n_all) args->gt[i] = 0;
									else args->gt[i] = (a1 << 4) | a2;
								}
							}
						}
						break;
						case BCF_BT_INT16:
						{
							if(fmt->n == 2) {
								int16_t *p = (int16_t *)(fmt->p + i * fmt->size);
								if(p[0] != bcf_int16_vector_end && p[1] != bcf_int16_vector_end) {
									int a1 = p[0] >> 1;
									int a2 = p[1] >> 1;
									if(a1 < 1 || a1 > n_all || a2 < 1 || a2 > n_all) args->gt[i] = 0;
									else args->gt[i] = (a1 << 4) | a2;
								}
							}
						}
						break;
						case BCF_BT_INT32:
						{
							if(fmt->n == 2) {
								int32_t *p = (int32_t *)(fmt->p + i * fmt->size);
								if(p[0] != bcf_int32_vector_end && p[1] != bcf_int32_vector_end) {
									int a1 = p[0] >> 1;
									int a2 = p[1] >> 1;
									if(a1 < 1 || a1 > n_all || a2 < 1 || a2 > n_all) args->gt[i] = 0;
									else args->gt[i] = (a1 << 4) | a2;
								}
							}
						}
						break;
						}
						if(args->gt[i]) passed = true;
					}
				}
			}
			if(passed) {
				kstring_t *s = ks_clear(&args->out_string);
				ksprintf(s, "%s\t%" PRId64 "\t%s", args->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + 1, id);
				for(int i = 0; i < ns; i++) {
					const int gt = args->gt[i];
					if(gt > 0) ksprintf(s, "\t%s%s", rec->d.allele[(gt >> 4) -1], rec->d.allele[(gt & 7) -1]);
					else kputs("\t00", s);
				}
				kputc('\n', s);
				htsFile * const fp = args->outfile;
				int r;
				if((fp)->format.compression != no_compression) r = bgzf_write((fp)->fp.bgzf, (s)->s, (s)->l);
				else r = hwrite((fp)->fp.hfile, (s)->s, (s)->l);
				if(r != (s)->l) error("output error writing to %s\n", (fp)->fn ? (fp)->fn : "<NULL>");
			}
		}
	}
	if(dbSNP_ctg) unload_dbSNP_ctg(dbSNP_ctg);
	bcf_destroy(rec);
	if(strcmp(args->outfilename, "-")) {
		hts_close(args->outfile);
		pthread_t md5_th;
		if(args->md5) pthread_create(&md5_th, NULL, md5_thread, args->outfilename);
		if(args->tabix) make_tabix_index(args->outfilename);
		if(args->md5) pthread_join(md5_th, NULL);
	}
}
