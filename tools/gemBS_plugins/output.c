#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "mextr.h"

int calc_phred(double z) {
  int phred;
  if(z <= 0.0) phred = 255;
  else {
    phred = (int)(-10.0 * log(z) / LOG10);
    if(phred > 255) phred = 255;
  }
  return phred;
}

static double *get_prob_dist(int ns, double *Q[]) {
  // Build up prob. distribution Q(i) where Q(i) = prob that i samples have genotype CG/CG
  double *p = Q[2];
  double *q0 = Q[0];
  double *q1 = Q[1];
  q0[0] = 1.0;
  for(int ix = 0; ix < ns; ix++) {
    double z = p[ix];
    q1[0] = q0[0] * (1.0 - z);
    for(int k = 1; k <= ix; k++) q1[k] = q0[k - 1] * z + q0[k] * (1.0 - z);
    q1[ix + 1] = q0[ix] * z;
    double *t = q0;
    q0 = q1;
    q1 = t;
  }
  return q0;
}

static char trans_base[256] = {
  ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A',
  ['Y'] = 'R', ['R'] = 'Y', ['S'] = 'S', ['W'] = 'W', ['K'] = 'M', ['M'] = 'K',
  ['B'] = 'V', ['V'] = 'B', ['D'] = 'H', ['H'] = 'D', ['N'] = 'N', ['.'] = '.'
};

void output_cpg(args_t *args, bcf1_t *rec, fmt_field_t *tags, gt_meth *sample_gt[], int idx, cpg_prob *cpg, double *Q[]) {
  static char *cx;
  static int32_t cx_n;
  static char *gt_iupac = "AMRWCSYGKT";
  static uint8_t gt_msk[] = {0x11, 0xb3, 0x55, 0x99, 0xa2, 0xf6, 0xaa, 0x54, 0xdc, 0x88};
  
  FILE *fp = args->cpgfile;
  int ns = bcf_hdr_nsamples(args->hdr);
  int min_n = args->min_num;
  int n1 = (int)(args->min_prop * (double)ns + 0.5);
  if(n1 > min_n) min_n = n1;
  if(fp != NULL) {
    // Build up prob. distribution Q(i) where Q(i) = prob that i samples have genotype CG/CG
    bool skip = true;
    for(int ix = 0; ix < ns; ix++) {
      gt_meth *g1 = sample_gt[idx]+ix, *g2 =sample_gt[idx ^ 1]+ix;
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
      int cx_len = bcf_get_info_values(args->hdr, rec, "CX", (void **)&cx, &cx_n, BCF_HT_STR);    
      int cx_sz = tags[FMT_CX].st[idx].ne / ns;
      if(args->mode == CPGMODE_COMBINED) {
	calc_cpg_meth(args, ns, cpg, sample_gt[idx], sample_gt[idx ^ 1]);
	fprintf(fp,"%s\t%d\t%d\t%.2s", args->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos, rec->pos + 2, cx_len >= 5 ? cx + 2 : ".");
	char *cx_p = tags[FMT_CX].st[idx].dat_p;
	int *mq_p1 = tags[FMT_MQ].st[idx].ne == ns ? tags[FMT_MQ].st[idx].dat_p : NULL; 
	int *mq_p2 = tags[FMT_MQ].st[idx^1].ne == ns ? tags[FMT_MQ].st[idx^1].dat_p : NULL;
	for(int ix = 0; ix < ns; ix++, cx_p += cx_sz) {
	  gt_meth *g1 = sample_gt[idx]+ix, *g2 =sample_gt[idx ^ 1]+ix;
	  if(!(g1->skip || g2->skip)) {
	    int gq = calc_phred(1.0 - exp(g1->gt_prob[g1->max_gt] + g2->gt_prob[g2->max_gt])); // Prob. of not being called genotype
	    fprintf(fp, "\t%c%c\tGQ=%d", gt_iupac[g1->max_gt], gt_iupac[g2->max_gt], gq);
	    if(g1->max_gt != 4 || g2->max_gt != 7) {
	      int dq = calc_phred(exp(g1->gt_prob[4] + g2->gt_prob[7])); // Prob. of being CG
	      fprintf(fp, ";DQ=%d", dq);
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
	    if(mq >= 0) fprintf(fp, ";MQ=%d", mq);
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
	    fprintf(fp, "\t%.3f\t%d\t%d\t%d\t%d", cpg[ix].m, ct[0], ct[1], ct[2], ct[3]);
	  } else {
	    fputs("\t.\t.\t.\t.\t.\t.\t.", fp);
	  }
	}
	fputc('\n', fp);
      } else {
	for(int pos = 0; pos < 2; pos++) {
	  int *mq_p = tags[FMT_MQ].st[idx ^ pos].ne == ns ? tags[FMT_MQ].st[idx ^ pos].dat_p : NULL; 
	  fprintf(fp,"%s\t%d\t%d\t%c", args->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + pos, rec->pos + pos + 1, cx_len >= 3 + pos ? cx[2 + pos] : '.');
	  char *cx_p = tags[FMT_CX].st[idx].dat_p;
	  for(int ix = 0; ix < ns; ix++, cx_p += cx_sz) {
	    gt_meth *g = sample_gt[idx ^ pos]+ix;
	    if(!g->skip) {
	      int gq = calc_phred(1.0 - exp(g->gt_prob[g->max_gt])); // Prob. of not being called genotype
	      fprintf(fp, "\t%c\tGQ=%d", gt_iupac[g->max_gt], gq);
	      if(g->max_gt != (pos ? 7 : 4)) {
		int dq = calc_phred(exp(g->gt_prob[pos ? 7 : 4])); // Prob. of being CG
		fprintf(fp, ";DQ=%d", dq);
	      }
	      int mq = -1;
	      if(mq_p != NULL) mq = mq_p[ix];
	      if(mq >= 0) fprintf(fp, ";MQ=%d", mq);
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
	      fprintf(fp, "\t%.3f\t%d\t%d\t%d\t%d", meth, ct[0], ct[1], ct[2], ct[3]);
	    } else {
	      fputs("\t.\t.\t.\t.\t.\t.\t.", fp);
	    }
	  }
	  fputc('\n', fp);
	}
      }
    }
  }
  if(args->output_noncpg) {
    fp=args->noncpgfile;
    assert(fp != NULL);
    for(int pos = 0; pos < 2; pos++) {
      for(int ix = 0; ix < ns; ix++) {
	double z = 0.0;
	gt_meth *g = sample_gt[idx ^ pos] + ix;
	if(!g->skip) {
	  if(!pos) {
	    if(g->counts[5] >= args->min_nc && (g->counts[5] + g->counts[7] >= args->min_inform)) {
	      if(args->sel_mode == SELECT_HOM) z = exp(g->gt_prob[4]);
	      else z = exp(g->gt_prob[1]) + exp(g->gt_prob[4]) + exp(g->gt_prob[5]) + exp(g->gt_prob[6]);
	      gt_meth *g2 = sample_gt[idx ^ 1] + ix;
	      z *= 1.0 - (exp(g2->gt_prob[2]) + exp(g2->gt_prob[5]) + exp(g2->gt_prob[7]) + exp(g2->gt_prob[8]));
	    }
	  } else {
	    if(g->counts[6] >= args->min_nc && (g->counts[6] + g->counts[4] >= args->min_inform)) {
	      if(args->sel_mode == SELECT_HOM) z = exp(g->gt_prob[7]);
	      else z = exp(g->gt_prob[2]) + exp(g->gt_prob[5]) + exp(g->gt_prob[7]) + exp(g->gt_prob[8]);
	      gt_meth *g2 = sample_gt[idx] + ix;
	      z *= 1.0 - (exp(g2->gt_prob[1]) + exp(g2->gt_prob[4]) + exp(g2->gt_prob[5]) + exp(g2->gt_prob[6]));
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
	int cx_len = bcf_get_info_values(args->hdr, rec, "CX", (void **)&cx, &cx_n, BCF_HT_STR);    
	int cx_sz = tags[FMT_CX].st[idx ^ pos].ne / ns;
	int *mq_p = tags[FMT_MQ].st[idx ^ pos].ne == ns ? tags[FMT_MQ].st[idx ^ pos].dat_p : NULL; 
	fprintf(fp,"%s\t%d\t%d\t%c", args->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + pos, rec->pos + pos + 1, cx_len >= 3 + pos ? cx[2 + pos] : '.');
	char *cx_p = tags[FMT_CX].st[idx ^ pos].dat_p;
	for(int ix = 0; ix < ns; ix++, cx_p += cx_sz) {
	  gt_meth *g = sample_gt[idx ^ pos] + ix;
	  if(!g->skip) {
	    int gq = calc_phred(1.0 - exp(g->gt_prob[g->max_gt])); // Prob. of not being called genotype
	    fprintf(fp, "\t%c\tGQ=%d", gt_iupac[g->max_gt], gq);
	    if(g->max_gt != (pos ? 7 : 4)) {
	      int dq = calc_phred(exp(g->gt_prob[pos ? 7 : 4])); // Prob. of being CG
	      fprintf(fp, ";DQ=%d", dq);
	    }
	    int mq = -1;
	    if(mq_p != NULL) mq = mq_p[ix];
	    if(mq >= 0) fprintf(fp, ";MQ=%d", mq);
	    if(!pos) {
	      if(cx_sz >= 5) fprintf(fp, ";CX=%.3s", cx_p + 2);
	    } else if(cx_sz >= 3) {
	      char tmp[3];
	      tmp[2] = trans_base[(int)cx_p[0]];
	      tmp[1] = trans_base[(int)cx_p[1]];
	      tmp[0] = trans_base[(int)cx_p[2]];
	      fprintf(fp, ";CX=%.3s", tmp);
	    }
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
	    fprintf(fp, "\t%g\t%d\t%d\t%d\t%d", meth, ct[0], ct[1], ct[2], ct[3]);
	  } else {
	    fputs("\t.\t.\t.\t.\t.\t.\t.\t.", fp);
	  }
	}
	fputc('\n', fp);
      }
    } 
  }
}							  

void output_nonconsec_noncpg(args_t *args, bcf1_t *rec, fmt_field_t *tags, gt_meth *sample_gt[], int idx, bool first, cpg_prob *cpg, double *Q[]) {
  static char *cx;
  static int32_t cx_n;
  static char *gt_iupac = "AMRWCSYGKT";
  static uint8_t gt_msk[] = {0x11, 0xb3, 0x55, 0x99, 0xa2, 0xf6, 0xaa, 0x54, 0xdc, 0x88};
  
  FILE *fp = args->cpgfile;
  int ns = bcf_hdr_nsamples(args->hdr);
  int min_n = args->min_num;
  int n1 = (int)(args->min_prop * (double)ns + 0.5);
  if(n1 > min_n) min_n = n1;
  fp=args->noncpgfile;
  assert(fp != NULL);
  for(int ix = 0; ix < ns; ix++) {
    double z = 0.0;
    gt_meth *g = sample_gt[idx] + ix;
    if(!g->skip) {
      if(first) {
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
    int cx_len = bcf_get_info_values(args->hdr, rec, "CX", (void **)&cx, &cx_n, BCF_HT_STR);    
    int cx_sz = tags[FMT_CX].st[idx].ne / ns;
    int *mq_p = tags[FMT_MQ].st[idx].ne == ns ? tags[FMT_MQ].st[idx].dat_p : NULL; 
    fprintf(fp,"%s\t%d\t%d\t%c", args->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos, rec->pos + 1, cx_len >= 3 ? cx[2] : '.');
    char *cx_p = tags[FMT_CX].st[idx].dat_p;
    for(int ix = 0; ix < ns; ix++, cx_p += cx_sz) {
      gt_meth *g = sample_gt[idx] + ix;
      if(!g->skip) {
	int gq = calc_phred(1.0 - exp(g->gt_prob[g->max_gt])); // Prob. of not being called genotype
	fprintf(fp, "\t%c\tGQ=%d", gt_iupac[g->max_gt], gq);
	if(g->max_gt != (first ? 4 : 7)) {
	  int dq = calc_phred(exp(g->gt_prob[first ? 4 : 7])); // Prob. of being CG
	  fprintf(fp, ";DQ=%d", dq);
	}
	int mq = -1;
	if(mq_p != NULL) mq = mq_p[ix];
	if(mq >= 0) fprintf(fp, ";MQ=%d", mq);
	if(cx_sz >= 5) fprintf(fp, ";CX=%.3s", cx_p + 2);
	int32_t ct[4];
	if(!first) {
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
	double meth = get_meth(g, !first);
	fprintf(fp, "\t%g\t%d\t%d\t%d\t%d", meth, ct[0], ct[1], ct[2], ct[3]);
      } else {
	fputs("\t.\t.\t.\t.\t.\t.\t.\t.", fp);
      }
    }
    fputc('\n', fp);
  }
}

static char *rgb_tab[11] = { "0,255,0", "55,255,0", "105,255,0", "155,255,0", "205,255,0", "255,255,0",
			     "255,205,0", "255,155,0", "255,105,0", "255,55,0", "255,0,0" };
  
void output_bedmethyl(args_t *args, bcf1_t *rec, fmt_field_t *tags, gt_meth *sample_gt[], int idx) {
  static char *cx;
  static int32_t cx_n,old_rid = 0xffffffff, old_pos = -1;
  
  if(rec->rid == old_rid && rec->pos <= old_pos) return;
  int ns = bcf_hdr_nsamples(args->hdr);
  if(ns > 1) return;
  gt_meth *g = sample_gt[idx];
  if(!g->skip) {
    int cx_len = bcf_get_info_values(args->hdr, rec, "CX", (void **)&cx, &cx_n, BCF_HT_STR);
    if(cx_len < 3) return;
    char strand;
    if(cx[2] == 'C') strand = '+';
    else if(cx[2] == 'G') strand = '-';
    else return;
    if(strand == '+' && cx_len < 5) return;
    char *cx_p = tags[FMT_CX].st[idx].dat_p;
    int cx_sz = tags[FMT_CX].st[idx].ne;
    char rtmp[8];
    if(strand == '+') {
      int k;
      for(k = 0; k < 3; k++) rtmp[k] = cx[k + 2];
      for(k = 0; k < 3 && k < cx_sz - 2; k++) rtmp[k + 4] = cx_p[k + 2];
      for(;k < 3; k++) rtmp[k + 4] = 'N';
    } else {
      int k;
      for(k = 0; k < 3; k++) rtmp[2 - k] = trans_base[(int)cx[k]];
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
    if(cov > 0 && args->wigfile) {
      FILE *fp = args->wigfile;
      if(fp != NULL) {
	if(rec->rid != old_rid) {
	  fprintf(fp, "variableStep chrom=%s\n", args->hdr->id[BCF_DT_CTG][rec->rid].key);
	}
	fprintf(fp, "%u\t%.4g\n", rec->pos + 1, 100.0 * m);
      }
    }
    FILE *fp = args->bedmethylfiles[btype];
    if(fp != NULL) {
      int gq = calc_phred(1.0 - exp(g->gt_prob[g->max_gt])); // Prob. of not being called genotype
      fprintf(fp, "%s\t%u\t%u\t\"%s\"\t%d\t%c\t%u\t%u\t%s\t%d\t%d\t%s\t%s\t%d\n",
	      args->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos, rec->pos + 1, args->bedmethyl_desc, cov > 1000 ? 1000 : cov, strand, 
	      rec->pos, rec->pos + 1, rgb_tab[(int)(m * 10.0 + 0.5)], cov, (int)(100.0 * m), rtmp, rtmp + 4, gq);
    }
    old_rid = rec->rid;
    old_pos = rec->pos;
  }
}
