/*
 * input_sam.c
 *
 *  Created on: Nov 29, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <sys/wait.h>
#include <pthread.h>
#include <htslib/sam.h>

#include "gem_tools.h"
#include "bs_call.h"

#define AA 0x0101
#define AC 0x0102
#define AG 0x0103
#define AT 0x0104
#define AN 0x0100
#define CA 0x0201
#define CC 0x0202
#define CG 0x0203
#define CT 0x0204
#define CN 0x0200
#define GA 0x0301
#define GC 0x0302
#define GG 0x0303
#define GT 0x0304
#define GN 0x0300
#define TA 0x0401
#define TC 0x0402
#define TG 0x0403
#define TT 0x0404
#define TN 0x0400
#define NA 0x0001
#define NC 0x0002
#define NG 0x0003
#define NT 0x0004
#define NN 0x0000

uint16_t bam_seq_decode[256] = {
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NA, AA, CA, NA, GA, NA, NA, NA, TA, NA, NA, NA, NA, NA, NA, NA,
	NC, AC, CC, NC, GC, NC, NC, NC, TC, NC, NC, NC, NC, NC, NC, NC,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NG, AG, CG, NG, GG, NG, NG, NG, TG, NG, NG, NG, NG, NG, NG, NG,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NT, AT, CT, NT, GT, NT, NT, NT, TT, NT, NT, NT, NT, NT, NT, NT,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN,
	NN, AN, CN, NN, GN, NN, NN, NN, TN, NN, NN, NN, NN, NN, NN, NN
};

static void get_seq_and_qual(align_details * const al, bam1_t *b) {
	bam1_core_t * const c = &b->core;
	int ix = (c->flag & BAM_FREVERSE) ? 1 : 0;
	uint32_t len = c->l_qseq;
	uint32_t len1 = 1 + ((len + 1) & ~(uint32_t)1);
	if (!al->read[ix]) {
		 al->read[ix] = gt_vector_new(len1, sizeof(uint8_t));
	} else {
		gt_vector_reserve(al->read[ix], len1, false);
	}
	const uint8_t *const sqb = bam_get_seq(b);
	uint8_t *const sq = gt_vector_get_mem(al->read[ix], uint8_t);
	uint16_t *tbuf = (uint16_t *)sq;
	for(int k = 0; k < (len + 1) >> 1; k++) tbuf[k] = bam_seq_decode[sqb[k]];
	const uint8_t *const ql = bam_get_qual(b);
	// Combine quality and base scores into an uint8_t;
	// Bits 0 and 1 code the bases A, C, G, T (0,1,2,3)
	// while bits 2-7 code the quality.  An N base is
	// set to have 0 quality, and by default is coded as 0x00.

	for(int k = 0; k < len; k++) {
		uint8_t q = ql[k];
		if(q > MAX_QUAL) q = MAX_QUAL;
		uint8_t c = sq[k];
		if(c) sq[k] = (c - 1) | (q << 2);
	}
	gt_vector_set_used(al->read[ix],len);
}

static uint32_t get_bam_misms(align_details * const al, bam1_t *b) {
	bam1_core_t * const c = &b->core;
	uint32_t reference_span = 0, position = 0;
	int ix = (c->flag & BAM_FREVERSE) ? 1 : 0;
	gt_vector_clear(al->mismatches[ix]);
	if(c->n_cigar) {
		const uint32_t * const cigar = bam_get_cigar(b);
		for(int i = 0; i < c->n_cigar; i++) {
			gt_misms misms;
			uint32_t length = bam_cigar_oplen(cigar[i]);
			switch(bam_cigar_opchr(cigar[i])) {
			case 'M':
			case '=':
			case 'X':
				position += length;
				reference_span += length;
				break;
		    case 'P': // Padding. Nothing specific implemented
		    case 'S': // Soft clipping. Looks similar to an insertion
		      misms.misms_type = SOFT;
		      misms.position = position;
		      misms.size = length;
		      position += length;
		      gt_vector_insert(al->mismatches[ix], misms, gt_misms);
		      break;
		    case 'H': // Hard clipping. Nothing specific implemented (we don't even store this)
		     break;
		    case 'I': // Insertion to the reference
		      misms.misms_type = DEL;
		      misms.position = position;
		      misms.size = length;
		      position += length;
		      gt_vector_insert(al->mismatches[ix], misms, gt_misms);
		      break;
		    case 'D': // Deletion from the reference
		      misms.misms_type = INS;
		      misms.position = position;
		      misms.size = length;
		      reference_span += length;
		      gt_vector_insert(al->mismatches[ix], misms, gt_misms);
		      break;
			}
		}
	}
	al->reference_span[ix]=reference_span;
	return position;
}

typedef enum { UNKNOWN_ALIGN, GEM, BOWTIE, NOVALIGN, BSMAP, BWAMETH } aligners;
static int aux_type2size[256] = {
		['A'] = 1, ['C'] = 1, ['c'] = 1, ['s'] = 2, ['S'] = 2, ['i'] = 4, ['I'] = 4,
		['f'] = 4, ['d'] = 8, ['Z'] = 'Z', ['H'] = 'H', ['B'] = 'B'
};

static gt_bs_strand get_bs_strand(bam1_t *b) {
	gt_bs_strand strand = NON_CONVERTED;
	uint8_t *s = bam_get_aux(b);
	bool ok = true;
	while(ok && (s + 4 <= b->data + b->l_data)) {
		 aligners align = UNKNOWN_ALIGN;
		 if(s[0] == 'Z') {
			 if(s[1] == 'B') align = NOVALIGN;
			 else if(s[1] == 'S') align = BSMAP;
		 } else if(s[0] == 'X') {
			 if(s[1] == 'G') align = BOWTIE;
			 else if(s[1] == 'B') align = GEM;
		 } else if(s[0] == 'Y' && s[1] == 'D') align = BWAMETH;
		 s+=2;
		 uint8_t type = *s++;
		 switch(type) {
		 case 'A':
			 if(align == GEM) {
				 if(*s == 'C') strand = STRAND_C2T;
				 else if(*s == 'G') strand = STRAND_G2A;
			 }
			 s++;
			 break;
		 case 'C':
		 case 'c':
			 s++;
			 break;
		 case 'S':
		 case 's':
			 if(s + 2 <= b->data + b->l_data) s+=2;
			 else ok = false;
			 break;
		 case 'I':
		 case 'i':
		 case 'f':
			 if(s + 4 <= b->data + b->l_data) s+=4;
			 else ok = false;
			 break;
		 case 'd':
			 if(s + 8 <= b->data + b->l_data) s+=8;
			 else ok = false;
			 break;
		 case 'Z':
			 if(align != UNKNOWN_ALIGN) {
				 if(align == BOWTIE || align == NOVALIGN) {
						if(*s == 'C') strand = STRAND_C2T;
						else if(*s == 'G') strand = STRAND_G2A;
				 } else if(align == BSMAP) {
						if(*s == '+') strand = STRAND_C2T;
						else if(*s == '-') strand = STRAND_G2A;
				 } else if(align == BWAMETH) {
						if(*s == 'f') strand = STRAND_C2T;
						else if(*s == 'r') strand = STRAND_G2A;
				 }
			 }
			 /* no break */
		 case 'H':
			 while (s<b->data + b->l_data && *s) s++;
			 if(s < b->data + b->l_data) s++;
			 else ok = false;
			 break;
		 case 'B': {
			 int sub_type_size = aux_type2size[(int)*(s++)];
			 if(s + 4 <= b->data + b->l_data && sub_type_size != 0) {
				 uint32_t n;
				 memcpy(&n, s, 4);
				 s+=4;
				 if(s + (n * sub_type_size) <= b->data + b->l_data) s += n * sub_type_size;
				 else ok = false;
			 } else ok = false;
		 }
		 break;
		 }
	}
	if(!ok) fprintf(stderr,"Error parsing SAM tags\n");
	return strand;
}

int get_next_align_details(htsFile * const sam_input, bam_hdr_t *hdr, hts_itr_t *itr, bam1_t *b, align_details * const al, const uint64_t thresh,
		const uint64_t max_template_len, bool keep_unmatched, bool ignore_dup, bool *reverse, gt_filter_reason * const filtered,
		uint32_t * const align_length, uint32_t * const alignment_flag) {

	*filtered = gt_flt_none;
	int ret;
	if(itr == NULL) ret = sam_read1(sam_input, hdr, b);
	else ret = sam_itr_next(sam_input, itr, b);
	if(ret >= 0) {
		ret = 0;
		bam1_core_t *c = &b->core;
		int flag = c->flag;
		if((flag & BAM_FPAIRED) && !keep_unmatched) {
			if (((flag & (BAM_FPROPER_PAIR | BAM_FUNMAP | BAM_FMUNMAP |
					BAM_FQCFAIL | BAM_FSECONDARY | BAM_FSUPPLEMENTARY |
					BAM_FDUP)) != BAM_FPROPER_PAIR)) {
				if(flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) *filtered = gt_flt_secondary;
				else if(flag & BAM_FUNMAP) *filtered = gt_flt_unmapped;
				else if(flag & BAM_FMUNMAP) *filtered = gt_flt_mate_unmapped;
				else if(flag & BAM_FQCFAIL) *filtered = gt_flt_qc;
				else if(flag & BAM_FDUP) {
					if(!ignore_dup) *filtered = gt_flt_duplicate;
				} else *filtered = gt_flt_not_correctly_aligned;
			}
		} else {
			if(flag & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FSECONDARY |
					BAM_FSUPPLEMENTARY | BAM_FDUP)) {
				if(flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) *filtered = gt_flt_secondary;
				else if(flag & BAM_FUNMAP) *filtered = gt_flt_unmapped;
				else if(flag & BAM_FQCFAIL) *filtered = gt_flt_qc;
				else if(flag & BAM_FDUP) *filtered = gt_flt_duplicate;
			}
		}
		bool mis_matched = (flag & (BAM_FMUNMAP | BAM_FPROPER_PAIR)) != BAM_FPROPER_PAIR;
		// Process flags
		*reverse = (flag & BAM_FREVERSE) ? true : false;
		const bool second_read = (flag & BAM_FREAD2);
		al->orientation = ((second_read && *reverse) || !(second_read || *reverse)) ? FORWARD : REVERSE;
		bool mult_seg = (flag & (BAM_FPAIRED | BAM_FMUNMAP)) == BAM_FPAIRED;
		if(*reverse) {
			al->forward_position = c->mpos + 1;
			al->reverse_position = c->pos + 1;
			al->mapq[1] = c->qual;
		} else {
			al->forward_position = c->pos + 1;
			al->reverse_position = c->mpos + 1;
			al->mapq[0] = c->qual;
		}
		if(c->qual < thresh && !(*filtered)) *filtered = gt_flt_mapq;
		*alignment_flag = flag;
		if(mult_seg) {
			if(c->tid != c->mtid) {
				if(!(*filtered)) *filtered = gt_flt_mismatch_chr;
				if(keep_unmatched) mis_matched = true;
			}
			if(!(*filtered)) {
				if(llabs(c->isize) > max_template_len) {
					*filtered = gt_flt_insert_size;
					if(keep_unmatched) mis_matched = true;
				}
			}
			if(*reverse) {
				if(c->pos < c->mpos) {
					if(!(*filtered)) *filtered = gt_flt_orientation;
					if(keep_unmatched) mis_matched = true;
				}
				if(mis_matched) al->forward_position = 0;
			} else {
				if(c->pos > c->mpos) {
					if(!(*filtered)) *filtered = gt_flt_orientation;
					if(keep_unmatched) mis_matched = true;
				}
				if(mis_matched) al->reverse_position = 0;
			}
		}
		 if(!mult_seg || mis_matched) *alignment_flag &= ~BAM_FPAIRED;
		 if(*filtered) {
			 if(!(keep_unmatched && (*filtered == gt_flt_insert_size || *filtered == gt_flt_mismatch_chr || *filtered == gt_flt_orientation))) ret = 1;
		 }
		 // We only decode the sequence and the tags if the read is not filtered
		 if(ret == 0) {
			 *align_length = get_bam_misms(al, b);
			 get_seq_and_qual(al, b);
			 al->bs_strand = get_bs_strand(b);
		 }

	} else {
		if(ret < -1) fprintf(stderr, "SAM input truncated\n");
	}
	return ret;
}
