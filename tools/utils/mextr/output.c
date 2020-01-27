#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>
#include <pthread.h>
#include <values.h>

#include "mextr.h"
#include "bbi.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/tbx.h"

void output_cpg(args_t *const args, rec_t ** const lrec, const int idx) {
	static char *gt_iupac = "AMRWCSYGKT";
	static uint8_t gt_msk[] = {0x11, 0xb3, 0x55, 0x99, 0xa2, 0xf6, 0xaa, 0x54, 0xdc, 0x88};
	double **Q = args->sample_Q;

	htsFile *fp = args->cpgfile;
	kstring_t *s = &args->pr_str[0];
	int ns = bcf_hdr_nsamples(args->hdr);
	int min_n = args->min_num;
	int n1 = (int)(args->min_prop * (double)ns + 0.5);
	if(n1 > min_n) min_n = n1;
	const rec_t * const rec1 = lrec[idx];
	const rec_t * const rec2 = lrec[idx ^ 1];
	if(fp != NULL) {
		// Build up prob. distribution Q(i) where Q(i) = prob that i samples have genotype CG/CG
		bool skip = true;
		for(int ix = 0; ix < ns; ix++) {
			gt_meth *g1 = rec1->sample_gt+ix, *g2 =rec2->sample_gt+ix;
			double z = 0.0;
			if(!(g1->skip || g2->skip)) {
				if((g1->counts[5] + g1->counts[7] >= args->min_inform) || (g2->counts[6] + g1->counts[4] >= args->min_inform)) {
					if(args->sel_mode == SELECT_HOM) {
						z = exp(g1->gt_prob[4] + g2->gt_prob[7]);
						if(g1->max_gt == 4 && g2->max_gt == 7) skip = false;
					} else {
						z = (exp(g1->gt_prob[1]) + exp(g1->gt_prob[4]) + exp(g1->gt_prob[5]) + exp(g1->gt_prob[6])) *
								(exp(g2->gt_prob[2]) + exp(g2->gt_prob[5]) + exp(g2->gt_prob[7]) + exp(g2->gt_prob[8]));
						if((g1->max_gt == 1 || (g1->max_gt >= 4 && g1->max_gt <= 6)) &&
								(g2->max_gt == 2 || g2->max_gt == 5 || g2->max_gt == 7 || g2->max_gt == 8)) skip = false;
					}
				}
			}
			Q[2][ix] = z;
		}
		double *p = get_prob_dist(ns, Q);
		double z = p[0];
		for(int i = 1; i <= ns && i < min_n; i++) z += p[i];
		int phred = calc_phred(z);
		if(!skip && phred >= args->sel_thresh) {
			int cx_sz = rec1->tags[FMT_CX].ne / ns;
			if(args->mode == CPGMODE_COMBINED) {
				calc_cpg_meth(args, ns, args->sample_cpg, lrec[idx]->sample_gt, lrec[idx ^ 1]->sample_gt);
				char ref[2] = {rec1->ref, rec2->ref};
				ksprintf(ks_clear(s), "%s\t%" PRId64 "\t%" PRId64 "\t%.2s", args->hdr->id[BCF_DT_CTG][rec1->rid].key, rec1->pos, rec1->pos + 2, ref);
				char *cx_p = rec1->tags[FMT_CX].dat_p;
				int *mq_p1 = rec1->tags[FMT_MQ].ne == ns ? rec1->tags[FMT_MQ].dat_p : NULL;
				int *mq_p2 = rec2->tags[FMT_MQ].ne == ns ? rec2->tags[FMT_MQ].dat_p : NULL;
				for(int ix = 0; ix < ns; ix++, cx_p += cx_sz) {
					gt_meth *g1 = rec1->sample_gt+ix, *g2 =rec2->sample_gt+ix;
					if(!(g1->skip || g2->skip)) {
						int gq = calc_phred(1.0 - exp(g1->gt_prob[g1->max_gt] + g2->gt_prob[g2->max_gt])); // Prob. of not being called genotype
						ksprintf(s, "\t%c%c\tGQ=%d", gt_iupac[g1->max_gt], gt_iupac[g2->max_gt], gq);
						if(g1->max_gt != 4 || g2->max_gt != 7) {
							int dq = calc_phred(exp(g1->gt_prob[4] + g2->gt_prob[7])); // Prob. of being CG
							ksprintf(s, ";DQ=%d", dq);
						}
						int mq = -1;
						if(mq_p1 != NULL) {
							if(mq_p2 != NULL) {
								double n1 = 0.0, n2 = 0.0;
								for(int k = 0; k < 8; k++) {
									n1 += (double)g1->counts[k];
									n2 += (double)g2->counts[k];
								}
								if(n1 + n2 > 0.0)  {
									double mq1 = (double)mq_p1[ix];
									double mq2 = (double)mq_p2[ix];
									mq = (int32_t)(0.5 + sqrt((mq1 * mq1 * n1 + mq2 * mq2 + n2) / (n1 + n2)));
								}
							} else mq = mq_p1[ix];
						} else if(mq_p2 != NULL) mq = mq_p2[ix];
						if(mq >= 0) ksprintf(s, ";MQ=%d", mq);
						int32_t ct[4];
						ct[0] = g1->counts[5] + g2->counts[6];
						ct[1] = g1->counts[7] + g2->counts[4];
						ct[2] = ct[3] = 0;
						uint8_t m = 1;
						uint8_t msk1 = gt_msk[g1->max_gt];
						uint8_t msk2 = gt_msk[g2->max_gt];
						for(int i = 0; i < 8; i++, m <<= 1) {
							ct[3] += g1->counts[i] + g2->counts[i];
							if(msk1 & m) ct[2] += g1->counts[i];
							if(msk2 & m) ct[2] += g2->counts[i];
						}
						ksprintf(s, "\t%.3f\t%d\t%d\t%d\t%d", args->sample_cpg[ix].m, ct[0], ct[1], ct[2], ct[3]);
					} else {
						kputs("\t.\t.\t.\t.\t.\t.\t.", s);
					}
				}
				kputc('\n', s);
			} else {
				for(int pos = 0; pos < 2; pos++) {
					rec_t *rec = lrec[idx ^ pos];
					int *mq_p = rec->tags[FMT_MQ].ne == ns ? rec->tags[FMT_MQ].dat_p : NULL;
					ksprintf(ks_clear(s), "%s\t%" PRId64 "\t%" PRId64 "\t%c", args->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + pos, rec->pos + pos + 1, rec->ref);
					char *cx_p = rec->tags[FMT_CX].dat_p;
					for(int ix = 0; ix < ns; ix++, cx_p += cx_sz) {
						gt_meth *g = rec->sample_gt+ix;
						if(!g->skip) {
							int gq = calc_phred(1.0 - exp(g->gt_prob[g->max_gt])); // Prob. of not being called genotype
							ksprintf(s, "\t%c\tGQ=%d", gt_iupac[g->max_gt], gq);
							if(g->max_gt != (pos ? 7 : 4)) {
								int dq = calc_phred(exp(g->gt_prob[pos ? 7 : 4])); // Prob. of being CG
								ksprintf(s, ";DQ=%d", dq);
							}
							int mq = -1;
							if(mq_p != NULL) mq = mq_p[ix];
							if(mq >= 0) ksprintf(s, ";MQ=%d", mq);
							int32_t ct[4];
							if(pos) {
								ct[0] = g->counts[6];
								ct[1] = g->counts[4];
							} else {
								ct[0] = g->counts[5];
								ct[1] = g->counts[7];
							}
							ct[2] = ct[3] = 0;
							uint8_t m = 1;
							uint8_t msk = gt_msk[g->max_gt];
							for(int i = 0; i < 8; i++, m <<= 1) {
								ct[3] += g->counts[i];
								if(msk & m) ct[2] += g->counts[i];
							}
							double meth = get_meth(g, pos);
							ksprintf(s, "\t%.3f\t%d\t%d\t%d\t%d", meth, ct[0], ct[1], ct[2], ct[3]);
						} else {
							kputs("\t.\t.\t.\t.\t.\t.\t.", s);
						}
					}
					kputc('\n', s);
				}
			}
			ks_output(fp, s);
		}
	}
}

void output_noncpg(args_t *const args, const rec_t * const rec) {
	if(!rec) return;
	static char *gt_iupac = "AMRWCSYGKT";
	static uint8_t gt_msk[] = {0x11, 0xb3, 0x55, 0x99, 0xa2, 0xf6, 0xaa, 0x54, 0xdc, 0x88};
	double **Q = args->sample_Q1;
	kstring_t *s = &args->pr_str[1];
	int ns = bcf_hdr_nsamples(args->hdr);
	int min_n = args->min_num;
	int n1 = (int)(args->min_prop * (double)ns + 0.5);
	if(n1 > min_n) min_n = n1;
	htsFile *fp = args->noncpgfile;
	assert(fp != NULL);
	bool cstrand = true;
	uint8_t gt = rec->max_common_gt;
	if(gt == 7 || gt == 2 || gt == 8) cstrand = false;
	else if(gt == 5 && rec->ref == 'G') cstrand = false;
	for(int ix = 0; ix < ns; ix++) {
		double z = 0.0;
		gt_meth *g = rec->sample_gt + ix;
		if(!g->skip) {
			if(cstrand) {
				if(g->counts[5] >= args->min_nc && (g->counts[5] + g->counts[7] >= args->min_inform)) {
					if(args->sel_mode == SELECT_HOM) z = exp(g->gt_prob[4]);
					else z = exp(g->gt_prob[1]) + exp(g->gt_prob[4]) + exp(g->gt_prob[5]) + exp(g->gt_prob[6]);
				}
			} else {
				if(g->counts[6] >= args->min_nc && (g->counts[6] + g->counts[4] >= args->min_inform)) {
					if(args->sel_mode == SELECT_HOM) z = exp(g->gt_prob[7]);
					else z = exp(g->gt_prob[2]) + exp(g->gt_prob[5]) + exp(g->gt_prob[7]) + exp(g->gt_prob[8]);
				}
			}
		}
		Q[2][ix] = z;
	}
	double *p = get_prob_dist(ns, Q);
	double z = p[0];
	for(int i = 1; i <= ns && i < min_n; i++) z += p[i];
	int phred = calc_phred(z);
	if(phred >= args->sel_thresh) {
		int cx_sz = rec->tags[FMT_CX].ne / ns;
		int *mq_p = rec->tags[FMT_MQ].ne == ns ? rec->tags[FMT_MQ].dat_p : NULL;
		ksprintf(ks_clear(s),"%s\t%"PRId64"\t%"PRId64"\t%c", args->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos, rec->pos + 1, rec->ref);
		char *cx_p = rec->tags[FMT_CX].dat_p;
		for(int ix = 0; ix < ns; ix++, cx_p += cx_sz) {
			gt_meth *g = rec->sample_gt + ix;
			if(!g->skip) {
				int gq = calc_phred(1.0 - exp(g->gt_prob[g->max_gt])); // Prob. of not being called genotype
				ksprintf(s, "\t%c\tGQ=%d", gt_iupac[g->max_gt], gq);
				if(g->max_gt != (cstrand ? 4 : 7)) {
					int dq = calc_phred(exp(g->gt_prob[cstrand ? 4 : 7])); // Prob. of being CG
					ksprintf(s, ";DQ=%d", dq);
				}
				int mq = -1;
				if(mq_p != NULL) mq = mq_p[ix];
				if(mq >= 0) ksprintf(s, ";MQ=%d", mq);
				if(cx_sz >= 5) ksprintf(s, ";CX=%.3s", cx_p + 2);
				int32_t ct[4];
				if(!cstrand) {
					ct[0] = g->counts[6];
					ct[1] = g->counts[4];
				} else {
					ct[0] = g->counts[5];
					ct[1] = g->counts[7];
				}
				ct[2] = ct[3] = 0;
				uint8_t m = 1;
				uint8_t msk = gt_msk[g->max_gt];
				for(int i = 0; i < 8; i++, m <<= 1) {
					ct[3] += g->counts[i];
					if(msk & m) ct[2] += g->counts[i];
				}
				double meth = get_meth(g, !cstrand);
				ksprintf(s, "\t%g\t%d\t%d\t%d\t%d", meth, ct[0], ct[1], ct[2], ct[3]);
			} else {
				kputs("\t.\t.\t.\t.\t.\t.\t.\t.", s);
			}
		}
		kputc('\n', s);
		ks_output(fp, s);
	}
}

static char *rgb_tab[11] = { "0,255,0", "55,255,0", "105,255,0", "155,255,0", "205,255,0", "255,255,0",
		"255,205,0", "255,155,0", "255,105,0", "255,55,0", "255,0,0" };

void output_bedmethyl(args_t *const args, const rec_t * const rec) {
	static int32_t old_rid = -1;
	static int64_t old_pos = -1;

	if(!rec) {
		// End of input
		// Flush buffers
		if(old_rid >= 0) {
			int id1 = args->id_trans[old_rid];
			assert(id1 >= 0);
			finish_bbi_blocks(args, id1);
		}
		return;
	}
	kstring_t *s = ks_clear(&args->pr_str[2]);
	if(rec->rid == old_rid && rec->pos <= old_pos) return;
	int ns = bcf_hdr_nsamples(args->hdr);
	if(ns > 1) return;
	gt_meth *g = rec->sample_gt;

	if(!g->skip) {
		char strand;
		if(rec->ref == 'C') strand = '+';
		else if(rec->ref == 'G') strand = '-';
		else return;

		char *cx_p = rec->tags[FMT_CX].dat_p;
		int cx_sz = rec->tags[FMT_CX].ne;
		char rtmp[8];
		if(strand == '+') {
			int k;
			for(k = 0; k < 3; k++) rtmp[k] = rec->cx[k + 2];
			for(k = 0; k < 3 && k < cx_sz - 2; k++) rtmp[k + 4] = cx_p[k + 2];
			for(;k < 3; k++) rtmp[k + 4] = 'N';
		} else {
			int k;
			for(k = 0; k < 3; k++) rtmp[2 - k] = trans_base[(int)rec->cx[k]];
			for(k = 0; k < 3 && k < cx_sz; k++) rtmp[6 - k] = trans_base[(int)cx_p[k]];
			for(;k < 3; k++) rtmp[6 - k] = 'N';
		}
		bedmethyl_type btype = BEDMETHYL_NONE;
		assert(rtmp[0] == 'C');
		if(rtmp[1] == 'G') {
			btype = BEDMETHYL_CPG;
			rtmp[2] = rtmp[6] = 0;
		} else {
			btype = rtmp[2] == 'G' ? BEDMETHYL_CHG : BEDMETHYL_CHH;
			rtmp[3] = rtmp[7] = 0;
		}
		if(btype != BEDMETHYL_NONE) {
			int32_t ct[2];
			if(strand == '-') {
				ct[0] = g->counts[6];
				ct[1] = g->counts[4];
			} else {
				ct[0] = g->counts[5];
				ct[1] = g->counts[7];
			}
			int32_t cov = ct[0] + ct[1];
			double m = cov > 0 ? (double)ct[0] / (double)cov : 0.0;

			if(cov > 0) {
				htsFile *fp = args->bedmethylfiles[btype - 1];
				assert(fp);
				int gq = calc_phred(1.0 - exp(g->gt_prob[g->max_gt])); // Prob. of not being called genotype
				ksprintf(s, "%s\t%"PRId64"\t%"PRId64"\t", args->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos, rec->pos + 1);
				size_t l = s->l;
				ksprintf(s, "\"%s\"\t%d\t%c\t%"PRId64"\t%"PRId64"\t%s\t%d\t%d\t%s\t%s\t%d\n", args->bedmethyl_desc, cov > 1000 ? 1000 : cov,
						strand, rec->pos, rec->pos + 1, rgb_tab[(int)(m * 10.0 + 0.5)], cov, (int)(100.0 * m), rtmp, rtmp + 4, gq);
				size_t sz_fields = s->l - l - 1;
				ks_output(fp, s);
				int id = args->id_trans[rec->rid];
				assert(id >= 0);
				bbi_ctg_data_t * const cdata = args->ctg_data + id;
				if(rec->rid != old_rid) {
					if(old_rid >= 0) {
						int id1 = args->id_trans[old_rid];
						assert(id1 >= 0);
						finish_bbi_blocks(args, id1);
					}
					old_rid = rec->rid;
					for(int j = 0; j < (args->strand_specific ? 5 : 4); j++) {
						for(int i = 0; i < ZOOM_LEVELS; i++) args->bb_global[j].res_end[i] = 0;
					}
				}

				//
				// First handle bigBed data
				//
				bbi_data_t * bdata = args->ctg_data[id].bbi_data + btype - 1;

				// Check if this is a new block
				if(bdata->n_items == 0) {
					bdata->bbuf.offset = 0;
					bdata->bbuf.start = rec->pos;
					bdata->bbuf.end = rec->pos + 1;
				}

				// Add data to current block
				uint32_t dt[3] = {id, rec->pos, rec->pos + 1};
				kstring_t *s1 = args->bb_global[btype - 1].buffer;
				args->bb_global[btype - 1].n_rec++;
				const uint32_t * zoom_scales = args->bb_global[btype - 1].zoom_scales;
				kputsn_((char *)dt, sizeof(uint32_t) * 3, s1);
				kputsn(s->s + l, sz_fields, s1);
				s1->l++;
				bdata->bbuf.end = rec->pos + 1;
				// Compress and write out block if full
				if(++(bdata->n_items) >= ITEMS_PER_SLOT) finish_bb_block(args, id, btype - 1);

				// Collect count data for zoom levels

				uint32_t * res_end = args->bb_global[btype - 1].res_end;
				uint32_t * res_size = args->bb_global[btype - 1].res_size;
				for(int i = 0; i < ZOOM_LEVELS; i++) {
					if(rec->pos >= res_end[i]) {
						res_size[i]++;
						res_end[i] = rec->pos + zoom_scales[i];
					}
				}

				//
				// And now bigWig data
				//
				const int k = (!args->strand_specific || strand == '+') ? 3 : 4;
				bdata = args->ctg_data[id].bbi_data + k;
				s = args->bb_global[k].buffer;
				zoom_scales = args->bb_global[k].zoom_scales;
				bw_rec_t * const bw_rec = bdata->bw_rec + bdata->n_items;
				bw_rec->start = rec->pos;
				bw_rec->val = 100.0 * m;
				// Compress and write out block if full
				if(++(bdata->n_items) >= BW_ITEMS_PER_SLOT) finish_bw_block(args, id, k - 3);

				// Collect data for zoom levels
				res_end = args->bb_global[k].res_end;
				res_size = args->bb_global[k].res_size;
				for(int i = 0; i < ZOOM_LEVELS; i++) {
					if(rec->pos >= res_end[i]) {
						res_size[i]++;
						res_end[i] = rec->pos + zoom_scales[i];
					}
				}
				//
				// And now collect detailed Zoom data (for both bigBed and bigWig
				//
				zoom_dt_t * const zd = &cdata->zoom_data;
				int ix = rec->pos >> 1;
				uint8_t base_type = btype;
				if(strand == '-') base_type |= 4;
				if(ct[0] > 0) {
					base_type |= 8;
					hts_resize(float, zd->val_ix + 1, &zd->val_size, &zd->val, 0);
					zd->val[zd->val_ix++] = 100.0 * m;
				}
				if(!(rec->pos & 1)) base_type <<= 4;
				zd->base_type[ix] |= base_type;
				old_pos = rec->pos;
			}
		}
	}
}

int add_bb_zoom_data_item(const uint32_t ct, const int ctg, const int ix, const uint32_t start, const uint32_t end, args_t * const args) {
	int ret = 0;
	if(ct > 0) {
		const float fct = (float)ct;
		ret = 1;
		bbi_data_t * const bdata = args->ctg_data[ctg].bbi_data + ix;
		// Check if this is a new block
		if(bdata->n_items == 0) {
			bdata->bbuf.offset = 0;
			bdata->bbuf.start = start;
			bdata->bbuf.end = end;
		}
		// Add data to current block
		uint32_t dt[4] = {ctg, start, end, ct};
		float dt1[4] = {1.0, 1.0, fct, fct};
		kstring_t *s = args->bb_global[ix].buffer;
		kputsn_((char *)dt, sizeof(uint32_t) * 4, s);
		kputsn_((char *)dt1, sizeof(float) * 4, s);
		bdata->bbuf.end = end;
		// Compress and write out block if full
		if(++(bdata->n_items) >= ITEMS_PER_SLOT) finish_bb_block(args, ctg, ix);
	}
	return ret;
}

int add_bw_zoom_data_item(bw_zrec_t * const bwr, const int ctg, const int ix, uint32_t scale, args_t * const args) {
	int ret = 0;
	if(bwr->count > 0) {
		ret = 1;
		bbi_data_t * const bdata = args->ctg_data[ctg].bbi_data + ix + 3;
		uint32_t end = bwr->end_base;
		uint32_t start = end - scale;
		// Check if this is a new block
		if(bdata->n_items == 0) {
			bdata->bbuf.offset = 0;
			bdata->bbuf.start = start;
			bdata->bbuf.end = end;
		}
		// Add data to current block
		uint32_t dt[4] = {ctg, start, end, bwr->count};
		float dt1[4] = {bwr->min, bwr->max, bwr->x, bwr->xsq};
		kstring_t *s = args->bb_global[ix + 3].buffer;
		kputsn_((char *)dt, sizeof(uint32_t) * 4, s);
		kputsn_((char *)dt1, sizeof(float) * 4, s);
		bdata->bbuf.end = end;
		// Compress and write out block if full
		if(++(bdata->n_items) >= BW_ITEMS_PER_SLOT) finish_bb_block(args, ctg, ix + 3);
	}
	return ret;
}

void add_bb_zrec(const uint32_t ct, const int ctg, const int zoom_level, const int ix, const uint32_t end, args_t * const args) {
	assert(zoom_level > 0);
	if(ct > 0) {
		bbi_zblock_t * const bzb = args->ctg_data[ctg].bbi_data[ix].zblock + zoom_level - 1;
		hts_resize(bb_zrec_t, bzb->ix + 1, &bzb->size, &bzb->bb_rec, 0);
		bb_zrec_t * const zrec = bzb->bb_rec + (bzb->ix++);
		zrec->end_base = end;
		zrec->count = ct;
	}
}

void add_bw_zrec(bw_zrec_t * const bwr, const int ctg, const int zoom_level, const int ix, uint32_t scale, args_t * const args) {
	assert(zoom_level > 0);
	if(bwr->count > 0) {
		bbi_zblock_t * const bzb = args->ctg_data[ctg].bbi_data[ix + 3].zblock + zoom_level - 1;
		hts_resize(bw_zrec_t, bzb->ix + 1, &bzb->size, &bzb->bw_rec, 0);
		bw_zrec_t * const zrec = bzb->bw_rec + (bzb->ix++);
		memcpy(zrec, bwr, sizeof(bw_zrec_t));
	}
}

void process_zoom_rec(int ctg, uint32_t pos, uint8_t x, uint32_t *count, uint32_t ct[ZOOM_LEVELS][3], uint32_t end[ZOOM_LEVELS][3],
		bw_zrec_t bwr[ZOOM_LEVELS][2], float * const val, int * const val_ix, args_t * const args) {

	const int ix = (x & 3) - 1;
	const int ix1 = (args->strand_specific && ((x & 4)) ? 1 : 0);
	const uint32_t * scales = args->bb_global[ix].zoom_scales;
	const uint32_t * scales1 = args->bb_global[ix1 + 3].zoom_scales;
	args->bb_global[ix].total_bases++;
	args->bb_global[ix1 + 3].total_bases++;
	float m = 0.0;
	if(x & 8) {
		m = val[(*val_ix)++];
		double m1 = m;
		if(m1 < args->bb_global[ix1 + 3].min_x) args->bb_global[ix1 + 3].min_x = m1;
		if(m1 > args->bb_global[ix1 + 3].max_x) args->bb_global[ix1 + 3].max_x = m1;
		args->bb_global[ix1 + 3].sum_x += m1;
		args->bb_global[ix1 + 3].sum_xsq += m1 * m1;
	} else args->bb_global[ix1 + 3].min_x = 0.0;
	for(int k = 0; k < ZOOM_LEVELS; k++) {
		// bigBed zoom data
		if(pos >= end[k][ix]) {
			if(k == 0) count[ix] += add_bb_zoom_data_item(ct[0][ix], ctg, ix, end[0][ix] - scales[0], end[0][ix], args);
			else add_bb_zrec(ct[k][ix], ctg, k, ix, end[k][ix], args);
			ct[k][ix] = 1;
			end[k][ix] = pos + scales[k];
		} else ct[k][ix]++;
		// bigWig zoom data
		if(pos >= bwr[k][ix1].end_base) {
			if(k == 0) count[ix1 + 3] += add_bw_zoom_data_item(&bwr[k][ix1], ctg, ix1, scales1[0], args);
			else add_bw_zrec(&bwr[k][ix1], ctg, k, ix1, scales1[k], args);
			bwr[k][ix1].end_base = pos + scales1[k];
			bwr[k][ix1].count = 1;
			bwr[k][ix1].min = bwr[k][ix1].max = bwr[k][ix1].x = m;
			bwr[k][ix1].xsq = m * m;
		} else {
			bwr[k][ix1].count++;
			bwr[k][ix1].x += m;
			bwr[k][ix1].xsq += m * m;
			if(m < bwr[k][ix1].min) bwr[k][ix1].min = m;
			else if(m > bwr[k][ix1].max) bwr[k][ix1].max = m;
		}
	}
}

void make_tabix_index(char * const fname) {
	tbx_conf_t conf = tbx_conf_bed;
	conf.line_skip = 1;
	tbx_index_build(fname, 0, &conf);
}

void *md5_thread(void *p) {
	calc_file_md5(p);
	return NULL;
}

void *cpg_thread(void *p) {
	args_t * args = p;
	rec_buffer_t * const rec_buf = &args->rec_buf;
	rec_t *lrec[2] = {NULL, NULL};
	const struct timespec wt = {0, 5000};
	uint64_t curr_idx = 0;
	int ix = 0;
	pthread_mutex_lock(&rec_buf->mut);
	while(true) {
		while(curr_idx - rec_buf->first_index >= REC_BUF_SIZE){
			pthread_mutex_unlock(&rec_buf->mut);
			nanosleep(&wt, NULL);
			pthread_mutex_lock(&rec_buf->mut);
		}
		rec_t *rec = rec_buf->buf[curr_idx - rec_buf->first_index];
		while(!args->rec_finished && !(rec->tasks & REC_READY)) {
			pthread_mutex_unlock(&rec_buf->mut);
			nanosleep(&wt, NULL);
			pthread_mutex_lock(&rec_buf->mut);
		}
		if(!(rec->tasks & REC_READY) && args->rec_finished) break;
		rec_t * const rec1 = lrec[ix ^ 1];
		if((rec->tasks & (RJ_OUTPUT_CPG | REC_SKIP)) == RJ_OUTPUT_CPG) {
			lrec[ix] = rec;
			if(rec1) {
				// Are the entries consecutive ?
				if(rec1->rid == rec->rid && rec1->pos + 1 == rec->pos) {
					pthread_mutex_unlock(&rec_buf->mut);
					output_cpg(args, lrec, ix ^ 1);
					pthread_mutex_lock(&rec_buf->mut);
				}
				rec1->tasks &= ~RJ_OUTPUT_CPG;
			}
		} else {
			lrec[ix] = NULL;
			rec->tasks &= ~RJ_OUTPUT_CPG;
			if(rec1) rec1->tasks &= ~RJ_OUTPUT_CPG;
		}
		curr_idx++;
		ix ^= 1;
	}
	if(lrec[ix ^ 1]) lrec[ix ^ 1]->tasks &= ~RJ_OUTPUT_CPG;
	pthread_mutex_unlock(&rec_buf->mut);
	hts_close(args->cpgfile);
	args->cpgfile = NULL;
	pthread_t md5_th;
	if(args->calc_md5) pthread_create(&md5_th, NULL, md5_thread, args->cpgfilename);
	if(args->tabix) make_tabix_index(args->cpgfilename);
	if(args->calc_md5) pthread_join(md5_th, NULL);
	return NULL;
}

void *output_thread(void *p) {
	thr_info_t *th = p;
	args_t * args = th->args;
	rec_buffer_t * const rec_buf = &args->rec_buf;
	const struct timespec wt = {0, 5000};
	uint64_t curr_idx = 0;
	pthread_mutex_lock(&rec_buf->mut);
	while(true) {
		while(curr_idx - rec_buf->first_index >= REC_BUF_SIZE){
			pthread_mutex_unlock(&rec_buf->mut);
			nanosleep(&wt, NULL);
			pthread_mutex_lock(&rec_buf->mut);
		}
		rec_t *rec = rec_buf->buf[curr_idx - rec_buf->first_index];
		while(!args->rec_finished && !(rec->tasks & REC_READY)) {
			pthread_mutex_unlock(&rec_buf->mut);
			nanosleep(&wt, NULL);
			pthread_mutex_lock(&rec_buf->mut);
		}
		if(!(rec->tasks & REC_READY) && args->rec_finished) break;
		if((rec->tasks & (th->job_flag | REC_SKIP)) == th->job_flag) {
			pthread_mutex_unlock(&rec_buf->mut);
			th->output(args, rec);
			pthread_mutex_lock(&rec_buf->mut);
		}
		curr_idx++;
		rec->tasks &= ~th->job_flag;
	}
	pthread_mutex_unlock(&rec_buf->mut);
	// Indicate end of input for bbi output
	th->output(args, NULL);
	if(th->job_flag == RJ_OUTPUT_NONCPG) {
		hts_close(args->noncpgfile);
		args->noncpgfile = NULL;
		args->cpgfile = NULL;
		pthread_t md5_th;
		if(args->calc_md5) pthread_create(&md5_th, NULL, md5_thread, args->noncpgfilename);
		if(args->tabix) make_tabix_index(args->noncpgfilename);
		if(args->calc_md5) pthread_join(md5_th, NULL);
	} else {
		for(int i = 0; i < 3; i++) {
			hts_close(args->bedmethylfiles[i]);
			args->bedmethylfiles[i] = NULL;
		}
		if(args->calc_md5) {
			pthread_t md5_th[3];
			for(int i = 0; i < 3; i++) pthread_create(md5_th + i, NULL, md5_thread, args->bedmethylnames[i]);
			for(int i = 0; i < 3; i++) pthread_join(md5_th[i], NULL);
		}
	}
	return NULL;
}

void *handle_bedmethyl_thread(void *p) {
	thr_info_t *th = p;
	args_t * args = th->args;
	pthread_t write_th, *compress_th;

	//
	// Main data output
	//
	int nb = args->compress_threads;
	compress_th = malloc(nb * sizeof(pthread_t));
	for(int i = 0; i < nb; i++) pthread_create(compress_th + i, NULL, bbi_compress_thread, args);
	pthread_create(&write_th, NULL, bbi_write_thread, args);
	output_thread(p);
	// Signal compression and write threads that there is no more input
	args->cblock_buf.end_of_input = true;
	// Wake up any compress or writing threads that are waiting for new input
	pthread_cond_broadcast(&args->cblock_buf.cond[1]);
	pthread_cond_broadcast(&args->cblock_buf.cond[2]);

	// Wait for compression and writing threads to complete
	for(int i = 0; i < args->compress_threads; i++) pthread_join(compress_th[i], NULL);
	pthread_join(write_th, NULL);
	// Store current file position of each file
	for(int i = 0; i < 3; i++) finish_bb_data_file(args, i);
	for(int i = 0; i < (args->strand_specific ? 2 : 1); i++) finish_bw_data_file(args, i);

	//
	// Create main index
	//
	fprintf(stderr, "Creating main indices\n");
	pthread_t bbi_index_thr[5];
	bbi_thr_info_t bbi_tinfo[5];
	int njobs = args->strand_specific ? 5 : 4;
	for(int ix = 0; ix < njobs; ix++) {
		FILE * const fp = ix < 3 ? args->bigbedfiles[ix] : args->bigwigfiles[ix - 3];
		fseek(fp, args->bbi_hdr[ix < 3 ? 0 : 1]->fullDataOffset, SEEK_SET);
		bbi_write(fp, args->bb_global[ix].n_rec);
		fseek(fp, args->bb_global[ix].index_offset, SEEK_SET);
		bbi_tinfo[ix].args = args;
		bbi_tinfo[ix].ix = ix;
		bbi_tinfo[ix].nrec = args->bb_global[ix].n_rec;
		pthread_create(bbi_index_thr + ix, NULL, bbi_create_index, bbi_tinfo + ix);
	}
	for(int i = 0; i < njobs; i++) pthread_join(bbi_index_thr[i], NULL);

	//
	// Create Zoom data & indices
	//

	// Initialize summary block
	for(int ix = 0; ix < njobs; ix++) {
		args->bb_global[ix].total_bases = 0;
		args->bb_global[ix].min_x = DBL_MAX;
		args->bb_global[ix].max_x = -DBL_MAX;
		args->bb_global[ix].sum_x = 0.0;
		args->bb_global[ix].sum_xsq = 0.0;
	}
	const bool tty = isatty(STDERR_FILENO);
	const uint64_t total_len = args->cumul_len[args->sr->regions->nseqs - 1];
	fprintf(stderr, "Creating zoom levels (%d)\n", ZOOM_LEVELS);
	for(int i = 0; i < ZOOM_LEVELS; i++) {
		clear_cblocks(args);;
		for(int ix = 0; ix < njobs; ix++) {
			FILE *fp = ix < 3 ? args->bigbedfiles[ix] : args->bigwigfiles[ix - 3];
			args->bb_global[ix].zoom_data_offset[i] = ftell(fp);
			bbi_write(fp, args->bb_global[ix].res_size[i]);
		}
		for(int k = 0; k < nb; k++) pthread_create(compress_th + k, NULL, bbi_compress_thread, args);
		pthread_create(&write_th, NULL, bbi_write_thread, args);
		uint32_t count[5] = {0, 0, 0, 0, 0};
		int old_complete = 0;
		uint64_t ctot = 0;
		for(int ctg = 0; ctg < args->sr->regions->nseqs; ctg++) {
			if(i == 0) {
				uint32_t end[ZOOM_LEVELS][3];
				uint32_t ct[ZOOM_LEVELS][3];
				bw_zrec_t bwr[ZOOM_LEVELS][2];
				for(int j = 0; j < ZOOM_LEVELS; j++) {
					for(int k = 0; k < 3; k++) end[j][k] = ct[j][k] = 0;
					for(int k = 0; k < 2; k++) memset(&bwr[j][k], 0, sizeof(bw_zrec_t));
				}
				zoom_dt_t * const zd = &args->ctg_data[ctg].zoom_data;
				int val_ix = 0;
				float * const val = zd->val;
				const uint8_t * tb = zd->base_type;
				uint32_t j_limit = (zd->len + 1) >> 1;
				for(uint32_t j = 0; j < j_limit; j++) {
					for(; !(*tb) && j < j_limit; j++, tb++);
					if(j == j_limit) break;
					const int pos = j << 1;
					const uint8_t x = *tb++;
					if(x & 0x30) process_zoom_rec(ctg, pos, x >> 4, count, ct, end, bwr, val, &val_ix, args);
					if(x & 3) process_zoom_rec(ctg, pos + 1, x, count, ct, end, bwr, val, &val_ix, args);
					if(tty) {
						const uint64_t tot = pos + 2 + ctot;
						int complete = (int)((tot * 100) / total_len);
						if(complete > old_complete) {
							old_complete = complete;
							fprintf(stderr,"Zoom level %d: %d%% complete   \r", i + 1, complete);
						}
					}
				}
				assert(val_ix == zd->val_ix);
				for(int ix = 0; ix < 3; ix++) {
					count[ix] += add_bb_zoom_data_item(ct[0][ix], ctg, ix, end[0][ix] - args->bb_global[ix].zoom_scales[0], end[0][ix], args);
					finish_bb_block(args, ctg, ix);
				}
				for(int ix = 0; ix < (args->strand_specific ? 2 : 1); ix++) {
					count[ix + 3] += add_bw_zoom_data_item(&bwr[0][ix], ctg, ix, args->bb_global[ix + 3].zoom_scales[0], args);
					finish_bb_block(args, ctg, ix + 3);
				}
				for(int k = 1; k < ZOOM_LEVELS; k++) {
					for(int ix = 0; ix < 3; ix++) {
						add_bb_zrec(ct[k][ix], ctg, k, ix, end[k][ix], args);
					}
					for(int ix = 0; ix < (args->strand_specific ? 2 : 1); ix++) {
						add_bw_zrec(&bwr[k][ix], ctg, k, ix, args->bb_global[ix + 3].zoom_scales[k], args);
					}
				}
				free(zd->base_type);
				if(zd->val) free(zd->val);
				zd->base_type = NULL;
				zd->val = NULL;
				zd->val_ix = zd->val_size = 0;
			} else {
				// Handle the subsequent zoom levels

				// For bigBed files
				for(int ix = 0; ix < 3; ix++) {
					bbi_zblock_t * const bzb = args->ctg_data[ctg].bbi_data[ix].zblock + i - 1;
					const uint32_t scale = args->bb_global[ix].zoom_scales[i];
					bb_zrec_t * zrec = bzb->bb_rec;
					if(zrec) {
						for(int j = 0; j < bzb->ix; j++, zrec++) {
							count[ix] += add_bb_zoom_data_item(zrec->count, ctg, ix, zrec->end_base - scale, zrec->end_base, args);
						}
						finish_bb_block(args, ctg, ix);
						free(bzb->bb_rec);
						bzb->bb_rec = NULL;
						bzb->ix = bzb->size = 0;
					}
				}
				// and for bigWig files
				for(int ix = 0; ix < (args->strand_specific ? 2 : 1); ix++) {
					bbi_zblock_t * const bzb = args->ctg_data[ctg].bbi_data[ix + 3].zblock + i - 1;
					const uint32_t scale = args->bb_global[ix + 3].zoom_scales[i];
					bw_zrec_t * zrec = bzb->bw_rec;
					if(zrec) {
						for(int j = 0; j < bzb->ix; j++, zrec++) {
							count[ix + 3] += add_bw_zoom_data_item(zrec, ctg, ix, scale, args);
						}
						finish_bb_block(args, ctg, ix + 3);
						free(bzb->bw_rec);
						bzb->bw_rec = NULL;
						bzb->ix = bzb->size = 0;
					}
				}
			}
			ctot = args->cumul_len[ctg];
			if(tty) {
				const uint64_t tot = ctot;
				int complete = (int)((tot * 100) / total_len);
				if(complete > old_complete) {
					old_complete = complete;
					fprintf(stderr,"Zoom level %d: %d%% complete\r", i + 1, complete);
				}
			}
		}
		// Signal compression and write threads that there is no more input
		args->cblock_buf.end_of_input = true;
		pthread_cond_broadcast(&args->cblock_buf.cond[1]);
		pthread_cond_broadcast(&args->cblock_buf.cond[2]);
		// Wait for compression and writing threads to complete
		for(int j = 0; j < args->compress_threads; j++) pthread_join(compress_th[j], NULL);
		pthread_join(write_th, NULL);
		for(int j = 0; j < njobs; j++) {
			assert(args->bb_global[j].res_size[i] == count[j]);
		}
		pthread_t bbi_index_thr[5];
		bbi_thr_info_t bbi_tinfo[5];
		for(int ix = 0; ix < njobs; ix++) {
			args->bb_global[ix].zoom_index_offset[i] = ftell(ix < 3 ? args->bigbedfiles[ix] : args->bigwigfiles[ix - 3]);
			bbi_tinfo[ix].args = args;
			bbi_tinfo[ix].ix = ix;
			bbi_tinfo[ix].nrec = count[ix];
			pthread_create(bbi_index_thr + ix, NULL, bbi_create_index, bbi_tinfo + ix);
		}
		for(int ix = 0; ix < njobs; ix++) pthread_join(bbi_index_thr[ix], NULL);
	}
	destroy_cblocks(args);
	free(compress_th);

	// Write headers, summary etc with correct values
	for(int ix = 0; ix < njobs; ix++) {
		FILE * const fp = (ix < 3 ? args->bigbedfiles[ix] : args->bigwigfiles[ix - 3]);
		bbi_header_t * const h = args->bbi_hdr[ix < 3 ? 0 : 1];
		h->fullIndexOffset = args->bb_global[ix].index_offset;
		h->uncompressBufSize = args->bb_global[ix].max_buf_size;
		if(ix < 3) {
			double x = args->bb_global[ix].total_bases;
			args->bb_global[ix].min_x = args->bb_global[ix].max_x = 1.0;
			args->bb_global[ix].sum_x = args->bb_global[ix].sum_xsq = x;
		}
		write_bbi_header(fp, h, args->bb_global + ix);
		fseek(fp, 0, SEEK_END);
		bbi_write(fp, h->magic);
		fclose(fp);
	}
	fprintf(stderr,"Zoom level %d: 100%% complete\nCalculating md5 sums of output files\n", ZOOM_LEVELS);
	if(args->calc_md5) {
		pthread_t md5_th[5];
		for(int ix = 0; ix < njobs; ix++) {
			char * const name = ix > 3 ? args->bigbednames[ix] : args->bigwignames[ix - 3];
			pthread_create(md5_th + ix, NULL, md5_thread, name);
		}
		for(int ix = 0; ix < njobs; ix++) pthread_join(md5_th[ix], NULL);
	}
	return NULL;
}
