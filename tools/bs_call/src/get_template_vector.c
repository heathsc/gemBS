/*
 * get_template_vector.c
 *
 *  Created on: Nov 24, 2019
 *      Author: heath
 */


#include <stdio.h>
#include <sys/wait.h>
#include <pthread.h>
#include <time.h>
#include <htslib/sam.h>

#include "gem_tools.h"
#include "bs_call.h"

static void handle_end_of_block(align_hash ** align_hash_p, gt_vector *align_list, const uint32_t max_pos, const int curr_tid,
		gt_vector * const free_hash_list, work_t * const work) {
	// Remove unmatched reads from hash
	align_hash *thash = *align_hash_p;
	while (thash) {
		align_hash *thash1 = thash->hh.next;
		add_al_hash_to_free_list(thash, free_hash_list);
		HASH_DEL(*align_hash_p, thash);
		thash = thash1;
	}
	int ix = gt_vector_get_used(align_list);
	if (ix) {
		pthread_mutex_lock(&work->process_mutex);
		while(work->align_list_waiting) {
			pthread_cond_wait(&work->process_cond2, &work->process_mutex);
		}
		pthread_mutex_unlock(&work->process_mutex);
		work->y_waiting = max_pos;
		int k = work->tid2id[curr_tid];
		assert(k >= 0);
		work->ctg_waiting = work->contigs[k];
		work->align_list_waiting = align_list;
		work->process_end = true;
		pthread_mutex_lock(&work->process_mutex);
		pthread_cond_signal(&work->process_cond1);
		pthread_mutex_unlock(&work->process_mutex);
	}
}

// Expects a SAM/BAM file sorted on coordinates.  Selects mapped pairs with TQ and
// TQ >= thresh, and make up vector of overlapping templates
gt_status read_input(htsFile *sam_input, gt_vector * align_list,sr_param *param) {
	GT_VECTOR_CHECK(align_list);
	gt_vector_clear(align_list);
	bam_hdr_t * const hdr = param->work.sam_header;
	int curr_tid = -1, old_tid = -1;
	bool chr_skip = false;
	bool new_contig = false;
	uint32_t max_pos = 0; // Position of rightmost end of current pileup
	uint32_t start_pos = 0; // Position of leftmost end of current pileup
	uint32_t read_idx = 0, curr_pos = 0, start_idx = 0;

	gt_string *tag = gt_string_new(128);
	align_details *al = NULL;
	align_hash *hash_base = NULL;
	gt_vector *free_list = gt_vector_new(256, sizeof(align_details *));
	gt_vector *free_hash_list = gt_vector_new(256, sizeof(align_hash *));
	gt_vector *al_hash_list = gt_vector_new(256, sizeof(align_hash *));
	const int n_reg = param->work.n_regions;
	int reg_ix = 0;
	hts_itr_t *itr = NULL;
	if(n_reg > 0 && param->work.sam_idx) {
		region_t * const reg = param->work.regions + reg_ix++;
		itr = sam_itr_queryi(param->work.sam_idx, reg->ctg->bam_tid, reg->start - 1, reg->stop);
		fprintf(stderr,"Processing region %s:%u-%u\n", reg->ctg->name, reg->start, reg->stop);
		param->work.curr_region = reg;
	}

	bam1_t *b = bam_init1();
	// Cycle through input lines until next read does not overlap with current pileup
	for(;;) {
		// Parse SAM Alignment
		bool reverse;
		gt_filter_reason filtered;
		uint32_t align_length;
		uint32_t alignment_flag;
		if (al == NULL) al = get_new_align_details(free_list);
		else clear_align_details(al);
		int ret;
		while(true) {
			ret = get_next_align_details(sam_input, hdr, itr, b, al, param->mapq_thresh, param->max_template_len,
					param->keep_unmatched, param->ignore_duplicates, &reverse, &filtered, &align_length, &alignment_flag);
			if(ret >= 0) break;
			if(ret != -1) return GT_STATUS_FAIL;
			if(itr != NULL) hts_itr_destroy(itr);
			if(reg_ix < n_reg && param->work.sam_idx) {
				region_t * const reg = param->work.regions + reg_ix++;
				itr = sam_itr_queryi(param->work.sam_idx, reg->ctg->bam_tid, reg->start - 1, reg->stop);
				fprintf(stderr,"Processing region %s:%u-%u\n", reg->ctg->name, reg->start, reg->stop);
				param->work.curr_region = reg;
			} else { // End of input
				handle_end_of_block(&hash_base, align_list, max_pos, curr_tid, free_hash_list, &param->work);
				return GT_STATUS_OK;
			}
		}
		if(ret > 0) {
			if(param->work.stats != NULL) {
				param->work.stats->filter_cts[filtered]++;
				param->work.stats->filter_bases[filtered] += b->core.l_qseq;
			}
			continue;
		}
		gt_string_set_nstring_static(tag, bam_get_qname(b), b->core.l_qname);
		bool new_block = false;
		new_contig = false;
		if(curr_tid < 0|| curr_tid != b->core.tid) {
			// New contig - must also be the start of a new block
			new_contig = new_block = true;
			chr_skip = false;
			old_tid = curr_tid;
			curr_tid = b->core.tid;
			assert(curr_tid >= 0);
			int k = param->work.tid2id[curr_tid];
			if(k < 0) chr_skip = true;
			else param->work.contigs[k]->curr_reg = param->work.curr_region;
			fprintf(stderr, "Processing chromosome %s (%s)\n", hdr->target_name[curr_tid], chr_skip ? "SKIP" : "OK");
		}

		// Test if alignment should be inserted into list, or whether this is a backwards facing member of a pair
		// in which case we should already have the mate

		// If the alignment is not paired or we (for some reason) only have information on one member of the pair,
		// we always add the read
		bool insert = true;
		if(!(chr_skip || new_contig)) {
			if((alignment_flag & BAM_FPAIRED) && al->forward_position > 0 && al->reverse_position > 0) {
				if (al->forward_position == al->reverse_position) {
					insert = ! check_for_tag_in_hash(hash_base, tag);
				} else if (reverse) insert = al->forward_position > al->reverse_position;
				else insert = al->forward_position < al->reverse_position;
			}
			// If we have a new read pair, check to see if we are still overlapping the current contig
			if (insert == true) {
				if(start_pos > 0) {
					if(al->forward_position > 0) {
						if(al->forward_position > max_pos && (al->reverse_position > max_pos || al->reverse_position == 0)) {
							if(al->forward_position - max_pos > 1) new_block = true;
						}
					} else if(al->reverse_position > max_pos && al->reverse_position - max_pos > 1) new_block = true;
				}
			}
		}
		// We need to start a new block, so we must process any remaining reads from the old block
		if (new_block == true) {
			// Remove unmatched reads from hash
			align_hash *thash = hash_base;
			while (thash) {
				align_hash *thash1 = thash->hh.next;
				add_al_hash_to_free_list(thash, free_hash_list);
				HASH_DEL(hash_base, thash);
				thash = thash1;
			}
			assert(insert == true);
			read_idx = 0;
			start_idx = 0;
			curr_pos = 0;
			int ix = gt_vector_get_used(align_list);
			// Put old block in processing queue
			if (ix) {
				align_details **al_p = gt_vector_get_mem(align_list, align_details *);
				uint32_t xx = (*al_p)->forward_position;
				if(xx == 0) xx = (*al_p)->reverse_position;
				pthread_mutex_lock(&param->work.process_mutex);
				while(param->work.align_list_waiting) {
					pthread_cond_wait(&param->work.process_cond2, &param->work.process_mutex);
				}
				pthread_mutex_unlock(&param->work.process_mutex);

				// If we are starting a new contig, make sure we use the previous
				// contig name for the last block!
				int tid = new_contig ? old_tid : curr_tid;
				assert(tid >= 0);
				int ix = param->work.tid2id[tid];
				assert(ix >= 0);
				param->work.ctg_waiting = param->work.contigs[ix];
				param->work.y_waiting = max_pos;
				gt_vector *new_free_list = param->work.free_list_waiting;
				param->work.free_list_waiting = NULL;
				param->work.align_list_waiting = align_list;
				pthread_mutex_lock(&param->work.process_mutex);
				pthread_cond_signal(&param->work.process_cond1);
				pthread_mutex_unlock(&param->work.process_mutex);
				if(new_free_list != NULL) {
					uint32_t used = gt_vector_get_used(free_list);
					ix = gt_vector_get_used(new_free_list);
					gt_vector_reserve(free_list, ix + used, false);
					memcpy(gt_vector_get_mem(free_list, align_details *)+used,
							gt_vector_get_mem(new_free_list, align_details *),
							sizeof(align_details *) * ix);
					gt_vector_set_used(free_list, ix + used);
					align_list = new_free_list;
					gt_vector_clear(align_list);
				} else align_list =  gt_vector_new(32, sizeof(align_details *));
			}
			if(new_contig && old_tid >= 0) {
				old_tid = -1;
				new_contig = false;
			}
			max_pos = start_pos = 0;
		}
		if(chr_skip) continue;

		// Update start and end position of contig
		uint32_t st, ml;
		if(reverse) {
			ml = al->reverse_position + al->reference_span[1];
			st = al->reverse_position;
		} else {
			ml = al->forward_position + al->reference_span[0];
			st= al->forward_position;
		}
		if(ml > max_pos) max_pos = ml;
		if(start_pos == 0 || start_pos > st) start_pos = st;
		uint32_t x;
		// Handle new read
		if(alignment_flag & BAM_FPAIRED) {
			// Backwards facing read, so the pair should already be present
			if (!insert) {
				align_hash *thash;
				HASH_FIND(hh, hash_base, gt_string_get_string(tag), gt_string_get_length(tag), thash);
				if (thash) {
					HASH_DEL(hash_base, thash);
					int ix = reverse ? 1 : 0;
					gt_vector *ts = thash->al->read[ix];
					thash->al->read[ix] = al->read[ix];
					al->read[ix] = ts;
					thash->al->mapq[ix] = al->mapq[ix];
					thash->al->reference_span[ix] = al->reference_span[ix];
					gt_vector *tv = thash->al->mismatches[ix];
					thash->al->mismatches[ix] = al->mismatches[ix];
					al->mismatches[ix] = tv;
					assert(al->forward_position == thash->al->forward_position && al->reverse_position == thash->al->reverse_position);
					gt_vector_set_elm(al_hash_list, thash->ix, align_hash *, NULL);
					add_al_hash_to_free_list(thash, free_hash_list);
				} else {
				  if(param->work.stats) {
				    param->work.stats->filter_cts[14]++;
				    param->work.stats->filter_bases[14] += gt_vector_get_used(al->read[reverse ? 1 : 0]);
				  }
					bool al_skip = false;
					// If the pair is missing, it could be because it was a duplicate and has
					// already been deleted.  In this case the read pair should lie within the current contig
					// In this case we can skip this read
					if(!param->keep_duplicates) {
						x = reverse ? al->reverse_position : al->forward_position;
						if(x >= start_pos) al_skip = true;
					}
					if(!al_skip) {
						if(param->keep_unmatched) {
							if(al->forward_position > 0) x = al->forward_position + align_length;
							else x = al->reverse_position + align_length;
							if(x > max_pos) max_pos = x;
							gt_vector_reserve(align_list, read_idx + 1, false);
							gt_vector_reserve(al_hash_list, read_idx + 1, false);
							if (gt_vector_get_used(align_list) <= read_idx) {
								gt_vector_set_used(align_list, read_idx + 1);
								gt_vector_set_used(al_hash_list, read_idx + 1);
							}
							gt_vector_set_elm(align_list, read_idx, align_details *, al);
							gt_vector_set_elm(al_hash_list, read_idx, align_hash *, NULL);
							read_idx++;
							al = NULL;
						} else {
							fprintf(stderr, "Warning not found: " PRIgts " %" PRIu32 " %" PRIu32 " %c\n",
									PRIgts_content(tag), al->forward_position, al->reverse_position,
									al->orientation == FORWARD ? '+' : '-');
						}
					}
				}
			} else {
				// Here we have a forward facing pair, so we need to store end to be
				// matched up later
				// But first we check for duplicates
				bool al_skip = false;
				if(!param->keep_duplicates) {
					uint32_t pos = al->forward_position > 0 ? al->forward_position : al->reverse_position;
					if(pos == curr_pos) {
						align_details **al_p = gt_vector_get_mem(align_list, align_details *);
						align_hash **alh_p = gt_vector_get_mem(al_hash_list, align_hash *);
						for(uint32_t ix = start_idx; ix < read_idx; ix++) {
							align_details *al1 = al_p[ix];
							if(al->forward_position == al1->forward_position && al->reverse_position == al1->reverse_position && al->bs_strand == al1->bs_strand) {
								int maxq = 0, maxq1 = 0;
								int kn = 0, kn1 = 0;
								for(int ix1 = 0; ix1 < 2; ix1++) {
									if(al->read[ix1] && gt_vector_get_used(al->read[ix1]) > 0) { maxq += al->mapq[ix1]; kn++; }
									if(al1->read[ix1] && gt_vector_get_used(al1->read[ix1]) > 0) { maxq1 += al1->mapq[ix1]; kn1++; }
								}
								maxq /= kn;
								maxq1 /= kn1;
								// If new read is better than previously stored read, replace read
								if((maxq1 < maxq) || (maxq == maxq1 && get_al_qual(al1) < get_al_qual(al))) {
									align_hash *thash;
									HASH_FIND(hh, hash_base, gt_string_get_string(tag), gt_string_get_length(tag), thash);
									gt_cond_fatal_error(thash && alh_p[ix], PARSE_SAM_DUPLICATE_SEQUENCE_TAG,PRIgts_content(tag));
									if(!thash) thash = alh_p[ix];
									al_p[ix] = al;
									if(thash != NULL) {
										HASH_DEL(hash_base, thash);
										make_align_hash(thash, tag, alignment_flag, ix, al);
									} else thash = make_new_align_hash(free_hash_list, tag, alignment_flag, ix, al);
									HASH_ADD_KEYPTR(hh, hash_base, gt_string_get_string(thash->tag), gt_string_get_length(thash->tag), thash);
									al = al1;
								}
								if(param->work.stats != NULL) {
									uint32_t len1 = al->read[0] ? gt_vector_get_used(al->read[0]) : 0;
									uint32_t len2 = al->read[1] ? gt_vector_get_used(al->read[1]) : 0;
									bool paired = len1 && len2;
									param->work.stats->filter_cts[gt_flt_duplicate] += paired ? 2 : 1;
									param->work.stats->filter_bases[gt_flt_duplicate] += len1 + len2;
								}
								al_skip = true;
							}
						}
					} else {
						curr_pos = pos;
						start_idx = read_idx;
					}
				}
				if(!al_skip) {
					align_hash *thash;
					HASH_FIND(hh, hash_base, gt_string_get_string(tag), gt_string_get_length(tag), thash);
					gt_cond_fatal_error(thash != NULL, PARSE_SAM_DUPLICATE_SEQUENCE_TAG,PRIgts_content(tag));
					thash = make_new_align_hash(free_hash_list, tag, alignment_flag, read_idx, al);
					HASH_ADD_KEYPTR(hh, hash_base, gt_string_get_string(thash->tag), gt_string_get_length(thash->tag), thash);
					gt_vector_reserve(align_list, read_idx + 1, false);
					gt_vector_reserve(al_hash_list, read_idx + 1, false);
					if (gt_vector_get_used(align_list) <= read_idx) {
						gt_vector_set_used(align_list, read_idx + 1);
						gt_vector_set_used(al_hash_list, read_idx + 1);
					}
					gt_vector_set_elm(align_list, read_idx, align_details *, al);
					gt_vector_set_elm(al_hash_list, read_idx, align_hash *, thash);
					read_idx++;
					al = NULL;
				}
			}
		} else { // Single (non-paired) reads
			bool al_skip = false;
			if(!param->keep_duplicates) {
				uint32_t pos = al->forward_position > 0 ? al->forward_position : al->reverse_position;
				if(pos == curr_pos) {
					align_details **al_p = gt_vector_get_mem(align_list, align_details *);
					align_hash **alh_p = gt_vector_get_mem(al_hash_list, align_hash *);
					for(uint32_t ix = start_idx; ix < read_idx; ix++) {
						align_details *al1 = al_p[ix];
						align_hash *thash = alh_p[ix];
						if(al->forward_position == al1->forward_position && al->reverse_position == al1->reverse_position &&
								al->bs_strand == al1->bs_strand && (thash == NULL || ((thash->alignment_flag & 9) == 9 || ((thash->alignment_flag & 9) == 0)))) {
							if((al1->mapq[0] < al->mapq[0]) || (al1->mapq[0] == al->mapq[0] && get_al_qual(al1)< get_al_qual(al))) {
								al_p[ix] = al;
								al = al1;
							}
							if(param->work.stats != NULL) {
								param->work.stats->filter_cts[gt_flt_duplicate]++;
								param->work.stats->filter_bases[gt_flt_none] += gt_vector_get_used(al->read[reverse ? 1 : 0]);
							}
							al_skip = true;
						}
					}
				} else{
					curr_pos = pos;
					start_idx = read_idx;
				}
			}
			if(!al_skip) {
				gt_vector_reserve(align_list, read_idx + 1, false);
				gt_vector_reserve(al_hash_list, read_idx + 1, false);
				if (gt_vector_get_used(align_list) <= read_idx) {
					gt_vector_set_used(align_list, read_idx + 1);
					gt_vector_set_used(al_hash_list, read_idx + 1);
				}
				gt_vector_set_elm(align_list, read_idx, align_details *, al);
				gt_vector_set_elm(al_hash_list, read_idx, align_hash *, NULL);
				read_idx++;
				al = NULL;
			}
		}
	}
	bam_destroy1(b);
	return GT_STATUS_OK;
}


