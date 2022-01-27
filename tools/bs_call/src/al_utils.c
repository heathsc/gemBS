/*
 * al_utils.c
 *
 *  Created on: Nov 24, 2019
 *      Author: heath
 */

#include <stdio.h>

#include "gem_tools.h"
#include "bs_call.h"

static int al_details, al_hash;

void print_mem() {
	fprintf(stderr,"%d %d\n", al_details, al_hash);
}

uint32_t get_al_qual(align_details *al) {
	uint32_t qual = 0, n = 0;
	for (int k = 0; k < 2; k++) {
		if(al->read[k]) {
			uint64_t rl = gt_vector_get_used(al->read[k]);
			uint8_t *sq = gt_vector_get_mem(al->read[k], uint8_t);
			for(int j = 0; j < rl; j++) {
				uint8_t q = GET_QUAL(sq[k]);
				if(q != FLT_QUAL) {
					qual += q;
					n++;
				}
			}
		}
	}
	return n > 0 ? qual / n : 0;
}

align_details *get_new_align_details(gt_vector * const free_list) {
	align_details *al = NULL;
	if (gt_vector_get_used(free_list)) {
		al = *(gt_vector_get_last_elm(free_list, align_details *));
		gt_vector_dec_used(free_list);
		clear_align_details(al);
	} else {
		al_details++;
		al = gt_malloc(sizeof(align_details));
		al->read[0] = al->read[1] = NULL;
		al->mismatches[0] = gt_vector_new(8, sizeof(gt_misms));
		al->mismatches[1] = gt_vector_new(8, sizeof(gt_misms));
		al->mapq[0] = al->mapq[1] = 0;
		al->forward_position = al->reverse_position = 0;
	}
	return al;
}

void clear_align_details(align_details * const al) {
	if (al->read[0]) {
		gt_vector_clear(al->read[0]);
		gt_vector_clear(al->mismatches[0]);
	}
	if (al->read[1]) {
		gt_vector_clear(al->read[1]);
		gt_vector_clear(al->mismatches[1]);
	}
	gt_vector_clear(al->mismatches[0]);
	gt_vector_clear(al->mismatches[1]);
	al->mapq[0] = al->mapq[1] = 0;
	al->forward_position = al->reverse_position = 0;
}

void make_align_hash(align_hash * const ah, gt_string * const tag, const uint32_t alignment_flag, const uint32_t ix, align_details * const al) {
	gt_string_copy(ah->tag, tag);
	ah->alignment_flag = alignment_flag;
	ah->ix = ix;
	ah->al = al;
}

align_hash *make_new_align_hash(gt_vector * const free_hash_list, gt_string * const tag, const uint32_t alignment_flag, const uint32_t ix, align_details * const al) {
	align_hash *ah = NULL;
	if(gt_vector_get_used(free_hash_list)) {
		ah = *(gt_vector_get_last_elm(free_hash_list, align_hash *));
		gt_vector_dec_used(free_hash_list);
	} else {
		al_hash++;
		ah = gt_malloc(sizeof(align_hash));
		ah->tag = gt_string_new(32);
	}
	make_align_hash(ah, tag, alignment_flag, ix, al);
	return ah;
}

void add_al_hash_to_free_list(align_hash * const ah, gt_vector * const free_hash_list) {
	gt_vector_reserve_additional(free_hash_list, 1);
	*(gt_vector_get_free_elm(free_hash_list, align_hash *)) = ah;
	gt_vector_inc_used(free_hash_list);
}

bool check_for_tag_in_hash(const align_hash * const hash, gt_string * const tag) {
	align_hash *thash;
	HASH_FIND(hh, hash, gt_string_get_string(tag), gt_string_get_length(tag), thash);
	return(thash != NULL);
}

static void left_trim_read(gt_vector *rd, const uint32_t l) {
	if(l > 0) {
		if(l >= gt_vector_get_used(rd)) gt_vector_clear(rd);
		else {
			uint32_t new_len = gt_vector_get_used(rd) - l;
			uint8_t *p = gt_vector_get_mem(rd, uint8_t);
			memmove(p, p + l, new_len);
			gt_vector_set_used(rd, new_len);
		}
	}
}

static void right_trim_read(gt_vector *rd, const uint32_t l) {
	if(l > 0) {
		if(l >= gt_vector_get_used(rd)) gt_vector_clear(rd);
		else gt_vector_get_used(rd) -= l;
	}
}

void trim_soft_clips(align_details * const al, bs_stats * const stats, uint32_t * const trim_left, uint32_t * const trim_right) {
	for (int k = 0; k < 2; k++) {
		if(al->read[k] == NULL) continue;
		uint32_t rl = gt_vector_get_used(al->read[k]);
		if(rl == 0) continue;
		uint32_t num_misms = gt_vector_get_used(al->mismatches[k]);
		// Trim Soft clips if present
		int nclip = 0;
		uint32_t adj = 0;
		for (int z = 0; z < num_misms; z++) {
			gt_misms *misms = gt_vector_get_elm(al->mismatches[k], z, gt_misms);
			if (misms->misms_type == SOFT) {
				if (z && z != num_misms - 1) gt_fatal_error_msg("Error in CIGAR - Soft clip not at extremity of read\n");
				nclip++;
				if (!misms->position) {
					if (misms->size >= rl)
						gt_fatal_error_msg("Error in CIGAR - Illegal soft clip (%d %" PRIu32 " %" PRIu32 " %" PRIu32 ")\n", z,
								misms->position, misms->size, rl);
					adj = misms->size;
					if(stats) stats->base_filter[base_clip] += adj;
					left_trim_read(al->read[k], adj);
					trim_left[k] = adj;
				} else {
					if (misms->position + misms->size != rl)
						gt_fatal_error_msg("Error in CIGAR - Illegal soft clip (%d %" PRIu32 " %" PRIu32 " %" PRIu32 ")\n", z,
								misms->position, misms->size, rl);
					right_trim_read(al->read[k], misms->size);
					trim_right[k] = misms->size;
					if(stats) stats->base_filter[base_clip] += misms->size;
				}
			} else if (nclip) {
				misms->position -= adj;
				gt_vector_set_elm(al->mismatches[k], z - nclip, gt_misms, *misms);
			}
		}
		if (nclip) {
			num_misms -= nclip;
			gt_vector_set_used(al->mismatches[k], num_misms);
		}
	}
}

void handle_overlap(align_details * const al, bs_stats * const stats, uint32_t * const trim_left, uint32_t * const trim_right) {
	uint32_t rdl[2]={0,0};
	if(al->read[0]) rdl[0]=gt_vector_get_used(al->read[0]);
	if(al->read[1]) rdl[1]=gt_vector_get_used(al->read[1]);
	bool rev;
	int32_t overlap;
	if(rdl[0] > 0 && rdl[1] > 0) {
		if (al->forward_position <= al->reverse_position) {
			overlap = al->reference_span[0] - al->reverse_position + al->forward_position;
			rev = false;
		} else {
			overlap = al->reference_span[1] + al->reverse_position - al->forward_position;
			rev = true;
		}
		// Look for overlapping reads - keep only best part
		if (al->forward_position + al->reference_span[0] >=	al->reverse_position) {
			// Find the overlapping parts
			uint32_t *rspan;
			// Find best quality read (normally read 1)
			rspan = al->reference_span;
			int trim_read;
			if (rspan[0] > rspan[1]) trim_read = 1;
			else if (rspan[0] < rspan[1]) trim_read = 0;
			else {
				uint32_t tot[2] = {0, 0};
				for (int k = 0; k < 2; k++) {
					const uint8_t * const qual = gt_vector_get_mem(al->read[k], uint8_t);
					int n = 0;
					for (int i = 0; i < rdl[k]; i++) {
						uint8_t q = GET_QUAL(qual[i]);
						if(q != FLT_QUAL) {
							tot[k] += q;
							n++;
						}
					}
					tot[k] = n > 0 ? tot[k] / n : 0;
				}
				trim_read = tot[0] <= tot[1] ? 0 : 1;
			}
			// Adjust starting position if a left trim is used
			if (!(rev && trim_read) && (rev || trim_read)) {
				if (trim_read) al->reverse_position += overlap;
				else al->forward_position += overlap;
			}
			// Find overlap point for read
			uint32_t num_misms = gt_vector_get_used(al->mismatches[trim_read]);
			if (!num_misms) {
				if ((rev && trim_read) || !(rev || trim_read)) {
					right_trim_read(al->read[trim_read], overlap);
				} else {
					left_trim_read(al->read[trim_read], overlap);
				}
			} else {
				bool trimmed = false;
				if ((rev && trim_read) || !(rev || trim_read)) {
					uint32_t xx = al->reference_span[trim_read] - overlap;
					int64_t adj = 0;
					int z;
					for (z = 0; z < num_misms; z++) {
						gt_misms *misms = gt_vector_get_elm(al->mismatches[trim_read], z, gt_misms);
						if (misms->position + adj >= xx) {
							int64_t trim = rdl[trim_read] - xx + adj;
							right_trim_read(al->read[trim_read], trim);
							num_misms = z;
							gt_vector_set_used(al->mismatches[trim_read], num_misms);
							trimmed = true;
							break;
						}
						if (misms->misms_type == INS) {
							if (misms->position + adj + misms->size >= xx) {
								int64_t trim = rdl[trim_read] - misms->position;
								misms->size = xx - (misms->position + adj);
								right_trim_read(al->read[trim_read], trim);
								num_misms = z + 1;
								gt_vector_set_used(al->mismatches[trim_read], num_misms);
								trimmed = true;
							}
							adj += misms->size;
						} else if (misms->misms_type == DEL) adj -= misms->size;
											}
					if (trimmed == false) {
						right_trim_read(al->read[trim_read], overlap);
					}
				} else {
					uint32_t xx = overlap;
					int z;
					int64_t adj = 0;
					for (z = 0; z < num_misms; z++) {
						gt_misms *misms = gt_vector_get_elm(al->mismatches[trim_read], z, gt_misms);
						if (misms->position + adj >= xx) {
							uint32_t trim = overlap - adj;
							left_trim_read(al->read[trim_read], trim);
							trimmed = true;
							if (z) {
								for (int z1 = z; z1 < num_misms; z1++) {
									gt_misms *misms = gt_vector_get_elm(al->mismatches[trim_read], z1, gt_misms);
									misms->position -= trim;
									gt_misms *misms1 = gt_vector_get_elm(al->mismatches[trim_read], z1 - z, gt_misms);
									gt_vector_set_elm(al->mismatches[trim_read], z1 - z,gt_misms, *misms);
									gt_vector_set_elm(al->mismatches[trim_read], z1,gt_misms, *misms1);
								}
								num_misms -= z;
								gt_vector_set_used(al->mismatches[trim_read], num_misms);
								break;
							} else {
								for (int z1 = 0; z1 < num_misms; z1++) {
									gt_misms *misms = gt_vector_get_elm(al->mismatches[trim_read], z1, gt_misms);
									misms->position -= trim;
								}
							}
							break;
						}
						if (misms->misms_type == INS) {
							if (misms->position + adj + misms->size >= xx) {
								misms->size = misms->position + misms->size + adj - xx;
								uint32_t trim = misms->position;
								left_trim_read(al->read[trim_read], trim);
								trimmed = true;
								int z2 = misms->size ? z : z + 1;
								for (int z1 = z2; z1 < num_misms; z1++) {
									gt_misms *misms = gt_vector_get_elm(al->mismatches[trim_read], z1, gt_misms);
									misms->position -= trim;
									if (z2) {
										gt_misms *misms1 = gt_vector_get_elm(al->mismatches[trim_read], z1 - z2, gt_misms);
										gt_vector_set_elm(al->mismatches[trim_read], z1 - z2,gt_misms, *misms);
										gt_vector_set_elm(al->mismatches[trim_read], z1, gt_misms, *misms1);
									}
								}
								num_misms -= z2;
								gt_vector_set_used(al->mismatches[trim_read], num_misms);
								break;
							}
							adj += misms->size;
						} else if (misms->misms_type == DEL) adj -= misms->size;
					}
					if (trimmed == false) {
						left_trim_read(al->read[trim_read], overlap - adj);
						gt_vector_set_used(al->mismatches[trim_read], 0);
					}
				}
			}
			uint32_t rdl1[2];
			rdl1[0] = al->read[0] ? gt_vector_get_used(al->read[0]) : 0;
			rdl1[1] = al->read[1] ? gt_vector_get_used(al->read[1]) : 0;
			if(stats) stats->base_filter[base_overlap] += (rdl[0] - rdl1[0] + rdl[1] - rdl1[1]);
			if ((rev && trim_read) || !(rev || trim_read)) { // Right trim
				trim_right[trim_read] += (rdl[trim_read] - rdl1[trim_read]);
			} else { // Left trim
				trim_left[trim_read] += (rdl[trim_read] - rdl1[trim_read]);
			}
			rdl[0] = rdl1[0];
			rdl[1] = rdl1[1];
		}
	}
}
