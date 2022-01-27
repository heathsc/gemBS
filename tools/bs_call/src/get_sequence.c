/*
 * get_sequence.c
 *
 *  Created on: Dec 2, 2019
 *      Author: heath
 */


#include <stdio.h>

#include <htslib/faidx.h>
#include "gem_tools.h"
#include "bs_call.h"

//static const char dectab[8] = {'N', 'A', 'C', 'G', 'T', 'N', 'N', 'N'};
static const char dectab[8] = {0, 1, 2, 3, 4, 0, 0, 0};

#define FILL_DECODE_BUF(buf, x) {uint16_t c = (x); for(int k = 4; k >= 0; k--) { buf[k] = dectab[(int)(c & 7)]; c >>= 3; }}

bool get_sequence_string(ctg_t * contig, const uint32_t x, const uint32_t sz, ctg_t * prev_ctg, gt_string * const ref, sr_param * const par) {

	bool err = false;
	if(sz == 0) return true;
	if(prev_ctg != contig && prev_ctg != NULL) free_sequence(prev_ctg);
	if(contig->seq == NULL) {
		err = load_sequence(contig, par->work.seq_idx, par->report_file != NULL);
		if(err) return err;
	}
	char *rp = ref->buffer;
	gt_string_resize(ref, sz + 1);
	uint32_t len = sz;
	uint32_t x1;
	for(x1 = x; x1 < contig->start_pos && len > 0; x1++, len--) *rp++ = 0;
	if(len > 0) {
		char decode[5];
		uint32_t off = x1 - contig->start_pos;
		uint16_t *p = contig->seq + off / 5;
		FILL_DECODE_BUF(decode, *p++);
		uint32_t k = off % 5;
		while(len > 0 && x1 < contig->end_pos) {
			*rp++ = decode[k++];
			if(k == 5) {
				FILL_DECODE_BUF(decode, *p++);
				k=0;
			}
			x1++;
			len--;
		}
		while(len--) *rp++ = 0;
	}
	*rp++ = 0;
	ref->length = rp - ref->buffer;
//	fprintf(stderr, "ref from %u %u %lu %s\n", x, sz, ref->length, err ? "BAD" : "GOOD");
	return err;
}
