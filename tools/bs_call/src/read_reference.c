/*
 * read_reference.c
 *
 *  Created on: Nov 24, 2019
 *      Author: heath
 */

#include <stdio.h>

#include <htslib/bgzf.h>
#include <htslib/khash.h>
#include <htslib/faidx.h>

#include "gem_tools.h"
#include "bs_call.h"

// Internal structures from htslib
typedef struct {
	int id;
    uint32_t line_len, line_blen;
    uint64_t len;
    uint64_t seq_offset;
    uint64_t qual_offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

struct __faidx_t {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
    enum fai_format_options format;
};

void free_sequence(ctg_t * const contig) {
	if(contig->seq != NULL) {
		free(contig->seq);
		contig->seq = NULL;
		contig->start_pos = 0;
		contig->end_pos = 0;
	}
}

bool load_sequence(ctg_t * const contig, faidx_t *idx1, const bool calc_gc) {

	static const uint8_t btab[256] = {
			['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4,
			['a'] = 1, ['c'] = 2, ['g'] = 3, ['t'] = 4,
	};
	static const int gc_tab[4] = {0, 1, 1, 0};

	gt_vector *gc_bins = calc_gc ? gt_vector_new(16374, sizeof(uint8_t)) : NULL;
	fprintf(stderr, "Loading reference for %s\n", contig->name);
	// Find offset in file
	struct __faidx_t *idx = (struct __faidx_t *)idx1;
	khiter_t iter = kh_get(s, idx->hash, contig->name);
	if(iter == kh_end(idx->hash)) return true; // Sequence not found
	faidx1_t v = kh_value(idx->hash, iter);
	fprintf(stderr,"len = %"PRIu64"\n", v.len);
	int ret = bgzf_useek(idx->bgzf, v.seq_offset, SEEK_SET);
	if(ret < 0) return true; // Seek failed
	// We pack the reference in 3 bits per base, 5 bases to a uint16_t;
	int c;
	uint32_t k = 0;
	uint16_t x = 0;
	int ct[2] = {0, 0};
	while(k < v.len) {
		c = bgzf_getc(idx->bgzf);
		if(c < 0) break;
		if(isgraph(c)) {
			int b = btab[c];
			if(b) {
				x = b;
				ct[gc_tab[b - 1]]++;
				break;
			}
			k++;
		}
	}
	contig->start_pos = k + 1;
	uint32_t len = v.len - k;
	uint32_t sz = (len + 9) / 5;
	uint16_t *p = contig->seq = gt_malloc(sz * sizeof(uint16_t));
	k = 1;
	int ix = 1;
	while(k < len) {
		c = bgzf_getc(idx->bgzf);
		if(c < 0) break;
		if(isgraph(c)) {
			int b = btab[c];
			x = (x << 3) | b;
			k++;
			if(!(k % 5)) {
				*p++ = x;
				x = 0;
			}
			if(calc_gc) {
				if(b) ct[gc_tab[b - 1]]++;
				ix++;
				if(ix == 100) {
					gt_vector_insert(gc_bins, (ct[0] + ct[1] == 100 ? ct[1] : 255), uint8_t);
					ix = ct[0] = ct[1] = 0;
				}
			}
		}
	}
	if(x) {
		x <<= (3 * (5 - (k % 5)));
		*p++ = x;
	}
	bool err = false;
	if(v.len != k + contig->start_pos - 1) { // Did not read in expected amount
		free(contig->seq);
		contig->seq = NULL;
		err = true;
	} else {
		contig->end_pos = v.len;
		if(calc_gc) {
			gt_ctg_stats * const ctg_stats = contig->ctg_stats;
			ctg_stats->nbins = gt_vector_get_used(gc_bins);
			if(ctg_stats->nbins) {
				ctg_stats->gc = malloc((size_t)ctg_stats->nbins);
				memcpy(ctg_stats->gc, gc_bins->memory, (size_t)ctg_stats->nbins);
			}
		}
		fprintf(stderr,"Read in %u bases on %s from %u - %u, stored in %zu bytes\n", contig->end_pos - contig->start_pos + 1,
				contig->name, contig->start_pos, contig->end_pos, sz * sizeof(uint16_t));
	}
	if(gc_bins != NULL) gt_vector_delete(gc_bins);
	return err;
}

bool get_sequence_index(sr_param * const par) {
	bool error = false;
	faidx_t *idx = fai_load(par->name_reference_file);
	if(idx != NULL) {
		fprintf(stderr, "Sequence index read in successfully\n");
		par->work.seq_idx = idx;
	} else error = true;
	return error;
}

