/*
 * print_vcf.c
 *
 *  Created on: Nov 24, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <unistd.h>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/khash.h>

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

#include "gem_tools.h"
#include "bs_call.h"

#define lfact(x) lfact2(x,lfact_store)

static void add_flt_counts(gt_vector *v, int ct, bool var) {
	if(ct >= v->elements_allocated) gt_vector_reserve(v, ct + 1, true);
	if(ct >= v->used) v->used = ct + 1;
	fstats_cts *c = gt_vector_get_elm(v, ct, fstats_cts);
	c->cts[var ? 1 : 0]++;
}

static dbsnp_ctg_t *dbSNP_ctg;

void _print_vcf_entry(bcf1_t *bcf, ctg_t * const ctg, gt_meth *gtm, const char *rf_ctxt,
		const uint32_t x, char *gt_store, const sr_param * const par) {
	static const char *ref_alt[10][5] = {
			{"A", "", "A", "A", "A"},       // AA
			{"AC", "C", "A", "AC", "AC"}, // AC
			{"AG", "G", "AG", "A", "AG"}, // AG
			{"AT", "T", "AT", "AT", "A"}, // AT
			{"C", "C", "", "C", "C"},       // CC
			{"CG", "CG", "G", "C", "CG"}, // CG
			{"CT", "CT", "T", "CT", "C"}, // CT
			{"G", "G", "G", "", "G"},       // GG
			{"GT", "GT", "GT", "T", "G"}, // GT
			{"T", "T", "T", "T", ""}        // TT
	};
	static const stats_mut mut_type[10][5] = {
			{mut_no, mut_no, mut_CA, mut_GA, mut_TA},   // AA
			{mut_no, mut_AC, mut_CA, mut_no, mut_no},       // AC
			{mut_no, mut_AG, mut_no, mut_GA, mut_no},       // AG
			{mut_no, mut_AT, mut_no, mut_no, mut_TA},       // AT
			{mut_no, mut_AC, mut_no, mut_GC, mut_TC},   // CC
			{mut_no, mut_no, mut_CG, mut_GC, mut_no},       // CG
			{mut_no, mut_no, mut_CT, mut_no, mut_TC},       // CT
			{mut_no, mut_AG, mut_CG, mut_no, mut_TG},   // GG
			{mut_no, mut_no, mut_no, mut_GT, mut_TG},       // GT
			{mut_no, mut_AT, mut_CT, mut_GT, mut_no},   // TT
	};
	static char *cs_str[10] = {"NA", "+", "-", "NA", "+",
			"+-", "+", "-", "-",  "NA"};
	static const int all_idx[10][5][2] = {
			{{1, 0}, {0, 0}, {1, 0}, {1, 0}, {1, 0}}, // AA
			{{1, 2}, {2, 0}, {1, 0}, {1, 2}, {1, 2}}, // AC
			{{1, 3}, {3, 0}, {1, 3}, {1, 0}, {1, 3}}, // AG
			{{1, 4}, {4, 0}, {1, 4}, {1, 4}, {1, 0}}, // AT
			{{2, 0}, {2, 0}, {0, 0}, {2, 0}, {2, 0}}, // CC
			{{2, 3}, {2, 3}, {3, 0}, {2, 0}, {2, 3}}, // CG
			{{2, 4}, {2, 4}, {4, 0}, {2, 4}, {2, 0}}, // CT
			{{3, 0}, {3, 0}, {3, 0}, {0, 0}, {3, 0}}, // GG
			{{3, 4}, {3, 4}, {3, 4}, {4, 0}, {3, 0}}, // GT
			{{4, 0}, {4, 0}, {4, 0}, {4, 0}, {0, 0}}  // TT
	};

	static const uint8_t gt_int[10][5] = {
			{0x44, 0x22, 0x44, 0x44, 0x44}, // AA
			{0x48, 0x24, 0x24, 0x48, 0x48}, // AC
			{0x48, 0x24, 0x48, 0x24, 0x48}, // AG
			{0x48, 0x24, 0x48, 0x48, 0x24}, // AT
			{0x44, 0x44, 0x22, 0x44, 0x44}, // CC
			{0x48, 0x48, 0x24, 0x24, 0x48}, // CG
			{0x48, 0x48, 0x24, 0x48, 0x24}, // CT
			{0x44, 0x44, 0x44, 0x22, 0x44}, // GG
			{0x48, 0x48, 0x48, 0x24, 0x24}, // GT
			{0x44, 0x44, 0x44, 0x44, 0x22}  // TT
	};
	static const char gt_flag[10][5] = {
			{0, 1, 0, 0, 0}, // AA
			{0, 0, 0, 0, 0}, // AC
			{0, 0, 0, 0, 0}, // AG
			{0, 0, 0, 0, 0}, // AT
			{0, 0, 0, 0, 0}, // CC
			{0, 0, 0, 0, 0}, // CG
			{0, 0, 0, 0, 0}, // CT
			{0, 0, 0, 0, 0}, // GG
			{0, 0, 0, 0, 0}, // GT
			{0, 0, 0, 0, 1}, // TT
	};
	static char *pbase = "NACGT";
	static char *iupac = "NAMRWCSYGKT";
	static int cflag[] = {0, 1, 0, 0, 1, 1, 1, 0, 0, 0};
	static int gflag[] = {0, 0, 1, 0, 0, 1, 0, 1, 1, 0};
	static ctg_t *old_ctg = NULL;
	static uint32_t old_x;
	char rs[512];
	char prf_ctxt[5];
	static uint32_t prev_cpg_x;
	static bool prev_cpg_flt;
	const double * const logp = par->defs.logp;
	const double * const lfact_store = par->defs.lfact_store;
	const bool * const gt_het = par->defs.gt_het;

	if(x == 0) return;
	if(old_ctg != ctg) {
		// Free GC bin list for old contig
		if(old_ctg && old_ctg->ctg_stats && old_ctg->ctg_stats->gc) {
			free(old_ctg->ctg_stats->gc);
			old_ctg->ctg_stats->gc = NULL;
			old_ctg->ctg_stats->nbins = 0;
		}
		old_ctg = ctg;
	} else if(x <= old_x) return;
	old_x = x;
	ctg_t *contig = par->work.vcf_ctg;
	gt_ctg_stats *const ctg_stats = contig->ctg_stats;
	uint64_t *counts = gtm->counts;
	uint32_t dp = 0, d_inf = 0, dp1 = 0;
	for (int i = 0; i < 4; i++) dp1 += counts[i];
	for (int i = 4; i < 8; i++) d_inf += counts[i];
	dp = dp1 + d_inf;
	if (!dp) return;
	size_t rs_len = 0;
	rs[0] = 0;
	uint8_t rs_found = false;
	if(dbSNP_ctg != NULL) rs_found = dbSNP_lookup_name(par->work.dbSNP_hdr, dbSNP_ctg, rs, &rs_len, x);
	for(int i = 0; i < 5; i++) prf_ctxt[i] = pbase[(int)rf_ctxt[i]];
	char rfc = prf_ctxt[2];
	int rfix = (int)rf_ctxt[2];
	int gt = gt_store[2] - 1;
	// Skip homozygous reference if AA or TT
	bool skip = (!par->all_positions && !(rs_found & 2) && gt_flag[gt][rfix]);
	double z = gtm->gt_prob[gt];
	int phred;
	double z1 = exp(z * LOG10);
	if (z1 >= 1.0)
		phred = 255;
	else {
		phred = (int)(-10.0 * log(1.0 - z1) / LOG10);
		if (phred > 255) phred = 255;
	}
	const char *alt = ref_alt[gt][rfix];
	const stats_mut mut = mut_type[gt][rfix];
	const int fs = (int)(-gtm->fisher_strand * 10.0 + 0.5);
	const uint32_t qd = dp1 > 0 ? phred / dp1 : phred;
	uint32_t flt = 0;
	if(!skip) {
		if(ctg->curr_reg) {
			skip = (x < ctg->curr_reg->start || x > ctg->curr_reg->stop);
		} else skip = x > ctg->end_pos;
	}
	if(!skip) {
		bcf_clear(bcf);
		kstring_t *str = &bcf->shared;
		// RID
		bcf->rid = ctg->vcf_rid;
		// POS
		bcf->pos = x - 1;
		// ID
		if(!rs_found) bcf_enc_size(str, 0, BCF_BT_CHAR);
		else {
			bcf_enc_size(str, rs_len, BCF_BT_CHAR);
			kputsn(rs, rs_len, str);
		}
		// REF
		bcf_enc_vchar(str, 1, &rfc);
		bcf->n_allele = 1;
		bcf->rlen = 1;
		// ALT
		while(*alt) {
			bcf->n_allele++;
			bcf_enc_vchar(str, 1, alt);
			alt++;
		}
		// QUAL
		bcf->qual = (float)phred;

		// FILTER
		if (phred < 20) flt |= 1;
		if (qd < 2) flt |= 2;
		if (fs > 60) flt |= 4;
		if (gtm->mq < 40) flt |= 8;
		int fid = par->work.vcf_ids[VCF_FLT_PASS];
		if(!flt) {
			bool mac1 = false;
			switch(gt) {
			case 1: // AC
				mac1 = (counts[1] + counts[5] + counts[7] <= 1 || counts[0] + counts[4] <= 1);
				break;
			case 2: // AG
				mac1 = (counts[2] + counts[6] <= 1 || counts[0] <= 1);
				break;
			case 3: // AT
				mac1 = (counts[3] + counts[7] <= 1 || counts[0] + counts[4] <= 1);
				break;
			case 5: // CG
				mac1 = (counts[2] + counts[6] + counts[4] <= 1 || counts[1] + counts[5] + counts[7] <= 1);
				break;
			case 6: // CT
				mac1 = (counts[3] <= 1 || counts[1] + counts[5] <= 1);
				break;
			case 8: // GT
				mac1 = (counts[3] + counts[7] <= 1 || counts[2] + counts[6] + counts[4] <= 1);
				break;
			}
			if(mac1) {
				flt |= 128;
				fid = par->work.vcf_ids[VCF_FLT_MAC1];
			}
		} else fid = par->work.vcf_ids[VCF_FLT_FAIL];
		bcf_enc_vint(str, 1, &fid, -1);
		// INFO field
		bcf->n_info = 1;
		bcf_enc_int1(str, par->work.vcf_ids[VCF_INFO_CX]);
		bcf_enc_vchar(str, 5, prf_ctxt);
		bcf->n_sample = 1;
	}

	// Genotype fields
	char ctxt[5];
	for (int i = 0; i < 5; i++)
		ctxt[i] = iupac[(int)gt_store[i]];
	char *cpg = ".";
	if ((gt_store[2] == 5 && gt_store[3] == 8) ||
			(gt_store[2] == 8 && gt_store[1] == 5))
		cpg = "CG";
	else if (gt_store[2] == 5) {
		if (gt_store[3]) {
			if (gflag[(int)gt_store[3] - 1])
				cpg = "H";
			else
				cpg = "N";
		} else
			cpg = "?";
	} else if (gt_store[2] == 8) {
		if (gt_store[1]) {
			if (cflag[(int)gt_store[1] - 1])
				cpg = "H";
			else
				cpg = "N";
		} else
			cpg = "?";
	} else if (cflag[(int)gt_store[2] - 1]) {
		if (gt_store[3]) {
			if (gflag[(int)gt_store[3] - 1])
				cpg = "H";
			else
				cpg = "N";
		} else
			cpg = "?";
	} else if (gflag[(int)gt_store[2] - 1]) {
		if (gt_store[1]) {
			if (cflag[(int)gt_store[1] - 1])
				cpg = "H";
			else
				cpg = "N";
		} else
			cpg = ".";
	}
	if(!skip) {
		bcf->n_fmt = 11;
		// Handle sample fields
		kstring_t *str = &bcf->indiv;
		// GT
		int32_t x[8];
		uint8_t gg = gt_int[gt][rfix];
		x[0] = (gg >> 4);
		x[1] = (gg & 0xf);
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_GT]);
		bcf_enc_vint(str, 2, x, 2);
		// FT
		char fbuf[24];
		int flen = 0;
		if(flt & 15) {
			char *p = fbuf;
			int f_ix = 0;
			uint32_t flt1 = flt & 31;
			bool first = true;
			while(flt1) {
				if(flt1 & 1) {
					if(!first) *p++ = ';';
					const char *p1 = par->defs.flt_name[f_ix];
					while((*p++ = *p1++));
					first = false;
				}
				flt1 >>= 1;
				f_ix++;
			}
			*p = 0;
			flen = p - fbuf;
		} else {
			strcpy(fbuf, "PASS");
			flen = 4;
		}
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_FT]);
		bcf_enc_size(str, flen, BCF_BT_CHAR);
		kputsn_(fbuf, flen, str);
		// DP
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_DP]);
		bcf_enc_int1(str, dp1);
		// MQ
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_MQ]);
		bcf_enc_int1(str, gtm->mq);
		// GQ
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_GQ]);
		bcf_enc_int1(str, phred);
		// QD
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_QD]);
		bcf_enc_int1(str, qd);
		// GL
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_GL]);
		const int *aix = all_idx[gt][rfix];
		float gtl[6];
		if (rfix) {
			int j = rfix * (9 - rfix) / 2 + rfix - 5;
			z = gtm->gt_prob[j];
			if(z < -99.999) z = -99.999;
		} else z = -99.999;
		gtl[0] = z;
		int n_gt = 1;
		for (int i = 0; i < 2 && aix[i] > 0; i++) {
			int j;
			if (rfix) {
				if (rfix < aix[i]) j = rfix * (9 - rfix) / 2 + aix[i] - 5;
				else j = aix[i] * (9 - aix[i]) / 2 + rfix - 5;
				z = gtm->gt_prob[j];
				if(z < -99.999) z = -99.999;
				gtl[n_gt++] = z;
			}
			for (int k = 0; k < i; k++) {
				if (aix[k] < aix[i]) j = aix[k] * (9 - aix[k]) / 2 + aix[i] - 5;
				else j = aix[i] * (9 - aix[i]) / 2 + aix[k] - 5;
				z = gtm->gt_prob[j];
			}
			j = aix[i] * (9 - aix[i]) / 2 + aix[i] - 5;
			z = gtm->gt_prob[j];
			if(z < -99.999) z = -99.999;
			gtl[n_gt++] = z;
		}
		bcf_enc_vfloat(str, n_gt, gtl);
		// MC8
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_MC8]);
		for(int i = 0; i < 8; i++) x[i] = counts[i];
		bcf_enc_vint(str, 8, x, -1);
		// AMQ
		int k = 0;
		for(int i = 0; i < 8; i++) if(counts[i] > 0) x[k++] = gtm->qual[i];
		if(k > 0) {
			bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_AMQ]);
			bcf_enc_vint(str, k, x, -1);
			bcf->n_fmt++;
		}
		// CS
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_CS]);
		int l = strlen(cs_str[gt]);
		bcf_enc_size(str, l, BCF_BT_CHAR);
		kputsn_(cs_str[gt], l, str);
		// CG
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_CG]);
		bcf_enc_size(str, 1, BCF_BT_CHAR);
		kputc_((int)*cpg, str);
		// CX
		bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_CX]);
		bcf_enc_size(str, 5, BCF_BT_CHAR);
		kputsn_(ctxt, 5, str);
		// FS
		if(gt_het[gt]) {
			bcf_enc_int1(str, par->work.vcf_ids[VCF_FMT_FS]);
			bcf_enc_int1(str, fs);
			bcf->n_fmt++;
		}
		htsFile *fp = par->work.vcf_file;
		if(bcf_write(fp, par->work.vcf_hdr, bcf)) gt_fatal_error_msg("Failed to write vcf/bcf record");
	}
	bs_stats *stats = par->work.stats;
	if(stats != NULL) {
		bool snp = false;
		bool multi = false;
		gt_cov_stats *gcov;
		HASH_FIND(hh, stats->cov_stats, &dp, sizeof(uint32_t), gcov);
		if(gcov == NULL) {
			gcov = calloc((size_t)1, sizeof(gt_cov_stats));
			gcov->coverage = dp;
			HASH_ADD(hh, stats->cov_stats, coverage, sizeof(uint32_t), gcov);
		}
		gcov->all++;
		int bn = (x - contig->start_pos) / 100;
		if(bn >= 0 && bn < ctg_stats->nbins) {
			int gc = ctg_stats->gc[bn];
			if(gc <= 100) gcov->gc_pcent[gc]++;
		}
		if(!skip) {
			if(alt[0] != '.') { //Variant site
				if(alt[1] == ',') multi = true;
				else snp = true;
				if(snp) {
					stats->snps[stats_all]++;
					ctg_stats->snps[stats_all]++;
					if(!flt) {
						stats->snps[stats_passed]++;
						ctg_stats->snps[stats_passed]++;
					}
				} else {
					stats->multi[stats_all]++;
					ctg_stats->multi[stats_all]++;
					if(!flt) {
						stats->multi[stats_passed]++;
						ctg_stats->multi[stats_passed]++;
					}
				}
				stats->qual[variant_sites][phred]++;
				gcov->var++;
			}
			add_flt_counts(stats->qd_stats, qd, gt_het[gt]);
			add_flt_counts(stats->fs_stats, fs, gt_het[gt]);
			add_flt_counts(stats->mq_stats, gtm->mq, gt_het[gt]);
			stats->filter_counts[gt_het[gt] ? 1 : 0][flt & 31]++;
			stats->qual[all_sites][phred]++;
			if(rs_found) {
				stats->dbSNP_sites[stats_all]++;
				ctg_stats->dbSNP_sites[stats_all]++;
				if(snp || multi) {
					stats->dbSNP_var[stats_all]++;
					ctg_stats->dbSNP_var[stats_all]++;
				}
				if(!flt) {
					stats->dbSNP_sites[stats_passed]++;
					ctg_stats->dbSNP_sites[stats_passed]++;
					if(snp || multi) {
						stats->dbSNP_var[stats_passed]++;
						ctg_stats->dbSNP_var[stats_passed]++;
					}
				}
			}
			if(!strcmp(cpg, "CG")) {
				bool ref_cpg = false;
				bool cpg_ok = false;
				uint32_t a,b;
				if(!strcmp(cs_str[gt], "+")) {
					prev_cpg_x = x;
					prev_cpg_flt = (flt != 0);
					if(!strncmp(prf_ctxt+2, "CG", 2)) ref_cpg = true;
					a = counts[5];
					b = counts[7];
					cpg_ok = true;
				} else if(!strcmp(cs_str[gt], "-")) {
					if(!strncmp(prf_ctxt+1, "CG", 2)) ref_cpg = true;
					if(x - prev_cpg_x == 1) {
						if(ref_cpg) {
							stats->CpG_ref[stats_all]++;
							ctg_stats->CpG_ref[stats_all]++;
							if(!(prev_cpg_flt || flt)) {
								stats->CpG_ref[stats_passed]++;
								ctg_stats->CpG_ref[stats_passed]++;
							}
						} else {
							ctg_stats->CpG_nonref[stats_all]++;
							stats->CpG_nonref[stats_all]++;
							if(!(prev_cpg_flt || flt)) {
								stats->CpG_nonref[stats_passed]++;
								ctg_stats->CpG_nonref[stats_passed]++;
							}
						}
					}
					a = counts[6];
					b = counts[4];
					cpg_ok = true;
				}
				if(cpg_ok) {
					if(ref_cpg) {
						stats->qual[CpG_ref_sites][phred]++;
					} else {
						stats->qual[CpG_nonref_sites][phred]++;
					}
					gcov->CpG[ref_cpg ? 0 : 1]++;
					gt_cov_stats *gcov1;
					HASH_FIND(hh, stats->cov_stats, &d_inf, sizeof(uint32_t), gcov1);
					if(gcov1 == NULL) {
						gcov1 = calloc((size_t)1, sizeof(gt_cov_stats));
						gcov1->coverage = d_inf;
						HASH_ADD(hh, stats->cov_stats, coverage, sizeof(uint32_t), gcov1);
					}
					gcov1->CpG_inf[ref_cpg ? 0 : 1]++;
					if(a + b) {
						double meth[101];
						double konst = lfact(a + b + 1) - lfact(a) - lfact(b);
						double sum = 0.0;
						if(a) meth[0] = 0.0;
						else sum = meth[0] = exp(konst);
						if(b) meth[100] = 0.0;
						else sum = (meth[100] = exp(konst));
						double da = (double)a;
						double db = (double)b;
						for(int i = 1; i < 100; i++) {
							sum += (meth[i] = exp(konst + logp[i - 1] * da + logp[99 - i] * db));
						}
						for(int i = 0; i < 101; i++) {
							double z = meth[i] / sum;
							if(ref_cpg) {
								stats->CpG_ref_meth[stats_all][i] += z;
								if(!flt) stats->CpG_ref_meth[stats_passed][i] += z;
							} else {
								stats->CpG_nonref_meth[stats_all][i] += z;
								if(!flt) stats->CpG_nonref_meth[stats_passed][i] += z;
							}
						}
					}
				}
			}
			if(mut != mut_no) {
				stats->mut_counts[mut][stats_all]++;
				if(!flt) stats->mut_counts[mut][stats_passed]++;
				if(rs_found) { // dbSNP
					stats->dbSNP_mut_counts[mut][stats_all]++;
					if(!flt) stats->dbSNP_mut_counts[mut][stats_passed]++;
				}
			}
		}
	}
}

static char gt_store[5];
uint32_t store_x = 0;
static gt_meth gtm_store[5];
static ctg_t *curr_ctg;
static char rf_ctxt[8];

// Print the last 2 entries in gt_store
void flush_vcf_entries(bcf1_t * bcf, const sr_param * const par) {
	if (curr_ctg && store_x) {
		for (int i = 0; i < 2; i++) {
			memmove(gt_store, gt_store + 1, 4);
			memmove(gtm_store, gtm_store + 1, 4 * sizeof(gt_meth));
			memmove(rf_ctxt, rf_ctxt + 1, 6);
			if (gt_store[2]) _print_vcf_entry(bcf, curr_ctg, gtm_store + 2, rf_ctxt, store_x - 1 + i, gt_store, par);
		}
		store_x = 0;
	}
}

void print_vcf_entry(bcf1_t *bcf, ctg_t * const ctg, gt_meth *gtm, const char *rf,
		const uint32_t x, const uint32_t xstart, bool skip, sr_param * const par) {
	if (curr_ctg != ctg) curr_ctg = NULL;
	if (curr_ctg == NULL) {
		curr_ctg = ctg;
		const char * const ctgname = ctg->name;
		if(par->work.dbSNP_hdr != NULL) {
			if(dbSNP_ctg) unload_dbSNP_ctg(dbSNP_ctg);
			HASH_FIND(hh, par->work.dbSNP_hdr->dbSNP, ctgname, strlen(ctgname), dbSNP_ctg);
			if(dbSNP_ctg) {
				bool ret = load_dbSNP_ctg(par->work.dbSNP_hdr, dbSNP_ctg);
				if(!ret) dbSNP_ctg = NULL;
			}
		}
	}
	uint32_t l = x - store_x;
	if (l < 5) {
		memmove(gt_store, gt_store + l, 5 - l);
		memmove(gtm_store, gtm_store + l, (5 - l) * sizeof(gt_meth));
		for (uint32_t i = 4; i >= 5 - l; i--) {
			gt_store[i] = 0;
		}
	} else memset(gt_store, 0, 5);
	assert(x > store_x);
	store_x = x;
	memcpy(gtm_store + 4, gtm, sizeof(gt_meth));
	if (x - xstart >= 4)
		strncpy(rf_ctxt, rf + x - xstart - 4, 7);
	else {
		uint32_t l = x - xstart;
		for (uint32_t i = 0; i < 4 - l; i++)
			rf_ctxt[i] = 0;
		strncpy(rf_ctxt + 4 - l, rf, 3 + l);
	}
	if(skip) gt_store[4] = 0;
	else {
		double z = gtm->gt_prob[0];
		int gt = 0;
		for (int i = 1; i < 10; i++)
			if (gtm->gt_prob[i] > z) {
				z = gtm->gt_prob[i];
				gt = i;
			}
		gt_store[4] = gt + 1;
	}
	if (gt_store[2]) _print_vcf_entry(bcf, curr_ctg, gtm_store + 2, rf_ctxt, x - 2, gt_store, par);
}

static char *scan_hdr_keys(char *tp, const int n, const char *keys[], char *c_p[], int ln_p[]) {
	for(int k = 0; k < 5; k++) c_p[k] = NULL;
	while(*tp && *tp != '\n') {
		char **s_p = NULL;
		int *l_p;
		for(int k = 0; k < n; k++) {
			if(!strncmp(tp, keys[k], 2) && tp[2] == ':') {
				s_p = c_p + k;
				l_p = ln_p + k;
				break;
			}
		}
		if(s_p != NULL) {
			tp += 3;
			*s_p = tp;
			*l_p = 0;
			while(*tp && *tp != '\n' && *tp != '\t') {
				tp++;
				(*l_p)++;
			}
		} else while(*tp && *tp != '\n' && *tp != '\t') tp++;
		if(*tp == '\t') tp++;
	}
	return tp;
}

void print_vcf_header(sr_param * const param, bam_hdr_t * hdr) {
	bcf_hdr_t *bh = bcf_hdr_init("w");
	if(bh == NULL) {
		fprintf(stderr, "Fatal error - out of memory in print_vcf_header\n");
		exit(-1);
	}
	param->work.vcf_hdr = bh;

	char mode[4];
	mode[0] = 'w';
	int i = 1;
	if(param->out_file_type & FT_BCF) {
		mode[i++] = 'b';
		if(param->out_file_type == FT_BCF) mode[i++] = 'u';
	} else if(param->out_file_type & FT_GZ) mode[i++] = 'z';
	mode[i] = 0;
	htsFile *hout = hts_open(param->output_file ? param->output_file : "-", mode);
	if(!hout) {
		fprintf(stderr, "Couldn't open output file\n");
		exit(-1);
	}
	param->work.vcf_file = hout;
	if(param->num_threads[OUTPUT_THREADS] > 0) hts_set_threads(hout, param->num_threads[OUTPUT_THREADS]);
	bcf_hdr_printf(bh, "##fileformat=%s", bcf_hdr_get_version(bh));
	gt_string *buf = gt_string_new(256);
	if(!param->benchmark_mode) {
		time_t cl = time(0);
		struct tm *tt = localtime(&cl);
		bcf_hdr_printf(bh, "##fileDate(dd/mm/yyyy)=%02d/%02d/%04d", tt->tm_mday, tt->tm_mon + 1, tt->tm_year + 1900);
		bcf_hdr_printf(bh, "##source=bs_call_v%s,under_conversion=%g,over_conversion=%g,mapq_thresh=%d,bq_thresh=%d", BS_CALL_VERSION, param->under_conv, param->over_conv, param->mapq_thresh, param->min_qual);
		if(param->work.dbSNP_hdr != NULL && param->work.dbSNP_hdr->dbSNP_header != NULL) bcf_hdr_printf(bh, "##dbsnp=<%s>", param->work.dbSNP_hdr->dbSNP_header);
		// Scan header lines for @RG line (Read Groups)
		// Keep track of barcodes encountered so we only print one lne per barcode
		typedef struct {
			char *str;
			UT_hash_handle hh;
		} sdict;
		sdict *bc_dict = NULL;
		char *c_p[5];
		int ln_p[5];
		const char *rg_keys[3] = {"BC", "SM", "DS"};
		char *tp = hdr->text;
		while(tp && *tp) {
			if(!(strncmp(tp, "@RG\t", 4))) {
				tp = scan_hdr_keys(tp + 4, 3, rg_keys, c_p, ln_p);
				if(c_p[0] != NULL) {
					sdict *sd = NULL;
					HASH_FIND(hh, bc_dict, c_p[0], ln_p[0], sd);
					if(sd == NULL) {
						sd = gt_alloc(sdict);
						sd->str = c_p[0];
						HASH_ADD_KEYPTR(hh, bc_dict, sd->str, ln_p[0], sd);
						gt_sprintf(buf, "##bs_call_sample_info=<ID=\"%.*s\"", ln_p[0], c_p[0]);
						if(c_p[1] != NULL) gt_sprintf_append(buf,",SM=\"%.*s\"", ln_p[1], c_p[1]);
						if(c_p[2] != NULL) gt_sprintf_append(buf,",DS=\"%.*s\"", ln_p[2], c_p[2]);
						gt_string_append_string(buf, ">", 1);
						bcf_hdr_append(bh, buf->buffer);
					}
				}
			}
			tp = strchr(tp, '\n');
			if(tp) tp++;
		}
		if(bc_dict != NULL) {
			sdict *tp, *tp1;
			HASH_ITER(hh, bc_dict, tp, tp1) {
				HASH_DEL(bc_dict,tp);
				free(tp);
			}
		}
	}
	// And now pick up the @SQ lines (Sequences)
	const char *sq_keys[5] = {"SN", "LN", "AS", "M5", "SP"};
	char *tp = hdr->text;
	while(tp && *tp) {
		char *c_p[5];
		int ln_p[5];
		if(!(strncmp(tp, "@SQ\t", 4))) {
			tp = scan_hdr_keys(tp + 4, 5, sq_keys, c_p, ln_p);
			if(c_p[0] != NULL && c_p[1] != NULL) {
				if (param->work.n_contigs > 0) {
					char *ts = gt_malloc(ln_p[0] + 1);
					memcpy(ts, c_p[0], ln_p[0]);
					ts[ln_p[0]] = 0;
					int tid = bam_name2id(hdr, ts);
					free(ts);
					if(tid < 0 || param->work.tid2id[tid] < 0) continue;
				}
				gt_sprintf(buf, "##contig=<ID=%.*s,length=%.*s", ln_p[0], c_p[0], ln_p[1], c_p[1]);
				if(c_p[2] != NULL) gt_sprintf_append(buf, ",assembly=%.*s", ln_p[2], c_p[2]);
				if(c_p[3] != NULL) gt_sprintf_append(buf, ",md5=%.*s", ln_p[3], c_p[3]);
				if(c_p[4] != NULL) gt_sprintf_append(buf, ",sp=%.*s", ln_p[4], c_p[4]);
				gt_string_append_string(buf, ">", 1);
				bcf_hdr_append(bh, buf->buffer);
			}
		}
		tp = strchr(tp, '\n');
		if(tp) tp++;
	}
	gt_string_delete(buf);
	bcf_hdr_append(bh, "##INFO=<ID=CX,Number=1,Type=String,Description=\"5 base sequence context (from position -2 to +2 on the positive strand) determined from the reference\">");
	bcf_hdr_append(bh, "##FILTER=<ID=fail,Description=\"No sample passed filters\">");
	bcf_hdr_append(bh, "##FILTER=<ID=q20,Description=\"Genotype Quality below 20\">");
	bcf_hdr_append(bh, "##FILTER=<ID=qd2,Description=\"Quality By Depth below 2\">");
	bcf_hdr_append(bh, "##FILTER=<ID=fs60,Description=\"Fisher Strand above 60\">");
	bcf_hdr_append(bh, "##FILTER=<ID=mq40,Description=\"RMS Mapping Quality below 40\">");
	bcf_hdr_append(bh, "##FILTER=<ID=mac1,Description=\"Minor allele count <= 1\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Sample Genotype Filter\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred scaled conditional genotype quality\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (non converted reads only)\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"RMS Mapping Quality\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=QD,Number=1,Type=Integer,Description=\"Quality By Depth (Variant quality / read depth (non-converted reads only))\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=MC8,Number=8,Type=Integer,Description=\"Base counts: non-informative for methylation (ACGT) followed by informative for methylation (ACGT)\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=AMQ,Number=.,Type=Integer,Description=\"Average base quailty for where MC8 base count non-zero\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=CS,Number=1,Type=String,Description=\"Strand of Cytosine relative to reference sequence (+/-/+-/NA)\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=CG,Number=1,Type=String,Description=\"CpG Status (from genotype calls: Y/N/H/?)\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=CX,Number=1,Type=String,Description=\"5 base sequence context (from position -2 to +2 on the positive strand) determined from genotype call\">");
	bcf_hdr_append(bh, "##FORMAT=<ID=FS,Number=1,Type=Integer,Description=\"Phred scaled log p-value from Fishers exact test of strand bias\"");
	if(param->sample_name) bcf_hdr_add_sample(bh, param->sample_name);
	if(bcf_hdr_write(hout, bh)) gt_fatal_error_msg("Failed to write vcf/bcf header");
	// Lookup correspondences between contigs and BCF rid and
	for(int i = 0; i < bh->n[BCF_DT_CTG]; i++) {
		int tid = bam_name2id(hdr, bh->id[BCF_DT_CTG][i].key);
		assert(tid >= 0);
		int id = param->work.tid2id[tid];
		assert(id >= 0);
		param->work.contigs[id]->vcf_rid = i;
	}
	// Lookup filter/info/fmt IDs
	vdict_t *d = (vdict_t *)bh->dict[BCF_DT_ID];
	khint_t k;
	char *names[16] = { "PASS", "fail", "mac1", "CX" ,"GT", "FT", "GL", "GQ", "DP", "MQ", "QD", "MC8", "AMQ", "CS", "CG", "FS"};

	for(int i = 0; i < 16; i++) {
		k = kh_get(vdict, d, names[i]);
		// Shouldn't happen because we've just added it!
		if(k == kh_end(d)) {
			fprintf(stderr, "Internal error - could not find ID %s in VCF header\n", names[i]);
			exit(-1);
		}
		param->work.vcf_ids[i] = kh_val(d, k).id;
	}
}
