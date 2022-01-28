/*
 * meth_profile.c
 *
 *  Created on: Nov 27, 2019
 *      Author: heath
 */

#include <stdio.h>

#include "gem_tools.h"
#include "bs_call.h"

// Set bit 2 for CA, CC or CT, and set bit 3 for AG, GG and TG
static const uint8_t rtab[64] = {
		0, 0, 0, 0, 0, 0, 0, 0,   // NX
		0, 0, 0, 8, 0, 0, 0, 0,   // AX
		0, 4, 4, 0, 4, 0, 0, 0,   // CX
		0, 0, 0, 8, 0, 0, 0, 0,   // GX
		0, 0, 0, 8, 0, 0, 0, 0,   // TX
		0, 0, 0, 0, 0, 0, 0, 0,   // NX
		0, 0, 0, 0, 0, 0, 0, 0,   // NX
		0, 0, 0, 0, 0, 0, 0, 0,   // NX
};

// Collect counts of conversion, mutation and non-conversion events in non-CpG contexts.  For each position in the read, we look
// at C and G positions not in CpG context (according to the reference).  We count C->C, C->T, G->G and G->A events separately on the different strands.
// We collect 4 counts per read position:
//    a = C->C on the G2A strand or G->G on the C2T strand (or either from a non-converted library)
//    b = C->T on the G2A strand or G->A on the C2T strand (or either from a non-converted library)
//    c = C->C on the C2T strand or G->G on the G2A strand
//    d = C->T on the C2T strand or G->A on the G2A strand
//
// We can use  b / (a + b) as an estimate of the sequencing error rate + mutation rate
// Assuming that non-CpG sites are truely non-methylated, then d / (c + d) is an estimate of conversion + sequencing error + mutation.
// so by combing both we can get an estimate of conversion per location on the read.
//
// We use a FSM to avoid branching.
// state has the previous and current reference bases as a 6 bit number, the high 3 bits code the previous base and the low 3 the current base
// 0 = N, 1 = A, 2 = C, 3 = G, 4 = T
// rtab[state] is either 4, 8 or 0.  If 4 then we have a C not followed by a G or an N, and if 8 we have a G not preceded by a C or an N
// The read vector has the combined quality and base information as a uint8_t.  The top 6 bits code the quality and the low 2 bits code the base.
// 0 = A, 1 = C, 2 = G, 3 = T.  N's have the quality set to 0 (and the base information is irrelevant)
// btab[] is a lookup table indexed on strand (as bits 8 and 9) combined with the read info above.  The encode value (xx) is a 4 bit number
// where the low 2 bits code for the count (0 = a, 1 = b, 2 = c, 3 = d), bit 2 corresponds to a potential valid count for a C/T base,
// while bit 3 corresponds to a potential valid count for a G/A base.  rtab[state] & btab[strand + read] will have either bit 2 or bit 3 set
// if we have a valid count.

void meth_profile(const align_details * const al, const uint32_t x, gt_vector * const orig_pos[2], const int max_pos, sr_param * const par) {
	if(!par->work.stats) return;
	gt_vector *mprof = par->work.stats->meth_profile;
	const char *ref_st = gt_string_get_string(par->work.ref1);
	if(max_pos + 1 > mprof->used) {
		gt_vector_reserve(mprof, max_pos + 1, true);
		gt_vector_set_used(mprof, max_pos + 1);
	}
	meth_cts * const mc = gt_vector_get_mem(mprof, meth_cts) + 1;	for(int k = 0; k < 2; k++) {
		if(!al->read[k]) continue;
		uint32_t rl = gt_vector_get_used(al->read[k]);
		if(rl == 0) continue;
		uint8_t *sp = gt_vector_get_mem(al->read[k], uint8_t);
		uint32_t pos = k ? al->reverse_position : al->forward_position;
		const uint8_t *rf = (uint8_t *)(ref_st + pos - x);
		const uint8_t * const btab = par->work.flt_tab + 256 * al->bs_strand;
		const int *posx = gt_vector_get_mem(orig_pos[k], int);
		uint8_t state = (pos > x ? ((ref_st[pos - x - 1] << 3) | (*rf++)) : 0);
		uint8_t mask = rtab[(int)state];
		for(uint32_t j = 0; j < rl; j++) {
			uint8_t xx = btab[(int)(*sp++)];
			uint64_t * const cts = mc[*posx++].conv_cts;
			uint8_t mask1 = (xx & mask) >> 1;
			state = (pos >= x ? (((state << 3) | (*rf++)) & 63) : 0);
			mask = rtab[(int)state];
			cts[xx & 3] += (((xx & mask) | mask1) >> 2) & 1;
		}
	}
}

