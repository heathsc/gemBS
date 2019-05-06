#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "mextr.h"

typedef struct {
  double e, k, ln_k[3];
} qual_prob;

static qual_prob q_prob[MAX_QUAL + 1];

void fill_base_prob_table(void) {
  for (int q = 0; q <= MAX_QUAL; q++) {
    double e = exp(-.1 * (double)q * LOG10);
    if(e > .5) e = .5;
    double k = e / (3.0 - 4.0 * e);
    q_prob[q].e = e;
    q_prob[q].k = k;
    q_prob[q].ln_k[0] = log(k);
    q_prob[q].ln_k[1] = log(0.5 + k);
    q_prob[q].ln_k[2] = log(1.0 + k);
  }
}

static inline void get_Z(double x1, double x2, double k1, double k2, double l, double t, double *Z) {
  double lpt = l + t;
  double lmt = l - t;
  double d = (x1 + x2) * lmt;
  // w = 1, p = 1
  double sinm = (x1 * (lpt + 2.0 * k2) - x2 * (2.0 - lpt + 2.0 * k1)) / d;
  if(sinm < -1.0) sinm = -1.0;
  else if(sinm > 1.0) sinm = 1.0;
  Z[0] = 0.5 * (lmt * sinm + 2.0 - lpt);
  // w = 1, p = 1/2
  sinm = (x1 * (2.0 + lpt + 4.0 * k2) - x2 * (2.0 - lpt + 4.0 * k1)) / d;
  if(sinm < -1.0) sinm = -1.0;
  else if(sinm > 1.0) sinm = 1.0;
  Z[1] = 0.5 * (lmt * sinm + 2.0 - lpt);
  // w = 1/2, p = 1
  sinm = (x1 * (lpt + 4.0 * k2) - x2 * (2.0 - lpt + 4.0 * k1)) / d;
  if(sinm < -1.0) sinm = -1.0;
  else if(sinm > 1.0) sinm = 1.0;
  Z[2] = 0.5 * (lmt * sinm + 2.0 - lpt);
}

static void add_bias(double *ll, char rf, double ref_bias) {
  double lrb = log(ref_bias);
  double lrb1 = log(0.5 * (1.0 + ref_bias));
  memset(ll, 0, sizeof(double) * 10);
  switch (rf) {
  case 'A':
    ll[0] = lrb;
    ll[1] = ll[2] = ll[3] = lrb1;
    break;
  case 'C':
    ll[4] = lrb;
    ll[1] = ll[5] = ll[6] = lrb1;
    break;
  case 'G':
    ll[7] = lrb;
    ll[2] = ll[5] = ll[8] = lrb1;
    break;
  case 'T':
    ll[9] = lrb;
    ll[3] = ll[6] = ll[8] = lrb1;
    break;
  }
}

// This function is taken from genotype_model.c in bs_call
// As far as possible the two functions should be kept in sync
// (Yes, a shared library would make more sense...to do)
void calc_gt_prob(gt_meth *gt, args_t *args, char rf) {
  qual_prob qp[8];
  for(int i = 0; i < 8; i++)  qp[i] = q_prob[gt->aqual[i]];
  double l = 1.0 - args->under_conv;
  double t = args->over_conv;
  double n[8];
  for (int i = 0; i < 8; i++) n[i] = (double)gt->counts[i];
  double ll[10];
  double ref_bias = args->ref_bias;
  // Add in prior from reference
  add_bias(ll, rf, ref_bias);
  if (n[0]) {
    ll[0] += n[0] * qp[0].ln_k[2]; // AA
    double tz = n[0] * qp[0].ln_k[1]; 
    ll[1] += tz; // AC
    ll[2] += tz; // AG
    ll[3] += tz; // AT
    tz = n[0] * qp[0].ln_k[0];
    ll[4] += tz;  // CC
    ll[5] += tz;  // CG
    ll[6] += tz;  // CT
    ll[7] += tz;  // GG
    ll[8] += tz;  // GT
    ll[9] += tz; // TT
  }
  if (n[1]) {
    ll[4] += n[1] * qp[1].ln_k[2]; // CC
    double tz = n[1] * qp[1].ln_k[1];
    ll[1] += tz; // AC
    ll[5] += tz; // CG
    ll[6] += tz; // CT
    tz = n[1] * qp[1].ln_k[0];
    ll[0] += tz;  // AA
    ll[2] += tz;  // AG
    ll[3] += tz;  // AT
    ll[7] += tz;  // GG
    ll[8] += tz;  // GT
    ll[9] += tz; // TT
  }
  if (n[2]) {
    ll[7] += n[2] * qp[2].ln_k[2]; // GG
    double tz = n[2] * qp[2].ln_k[1];
    ll[2] += tz; // AG
    ll[5] += tz; // CG
    ll[8] += tz; // TG
    tz = n[2] * qp[2].ln_k[0];
    ll[0] += tz;  // AA
    ll[1] += tz;  // AC
    ll[3] += tz;  // AT
    ll[4] += tz;  // CC
    ll[6] += tz;  // CT
    ll[9] += tz; // TT
  }
  if (n[3]) {
    ll[9] += n[3] * qp[3].ln_k[2]; // TT
    double tz = n[3] * qp[3].ln_k[1];
    ll[3] += tz; // AT
    ll[6] += tz; // CT
    ll[8] += tz; // GT
    tz = n[3] * qp[3].ln_k[0];
    ll[0] += tz; // AA
    ll[1] += tz; // AC
    ll[2] += tz; // AG
    ll[4] += tz; // CC
    ll[5] += tz; // CG
    ll[7] += tz; // GG
  }
  double Z[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  if (n[5] + n[7] > 0.0) {
    get_Z(n[5], n[7], qp[5].k, qp[7].k, l, t, Z);
    for(int k = 0; k < 3; k++) gt->cmeth[k] = (Z[k] - 1.0 + l) / (l - t);
  }
  if (n[4] + n[6] > 0.0) {
    get_Z(n[6], n[4], qp[6].k, qp[4].k, l, t, Z+3);
    for(int k = 0; k < 3; k++) gt->gmeth[k] = (Z[k + 3] - 1.0 + l) / (l - t);
  }
  if (n[4]) {
    ll[0] += n[4] * qp[4].ln_k[2]; // AA
    ll[2] += log(1.0 - 0.5 * Z[4] + qp[4].k) * n[4]; // AG
    ll[7] += log(1.0 - Z[3] + qp[4].k) * n[4]; // GG
    double tz = log(0.5 * (1.0 - Z[5]) + qp[4].k) * n[4];
    ll[5] += tz; // CG
    ll[8] += tz; // GT
    tz = n[4] * qp[4].ln_k[1];
    ll[1] += tz; // AC
    ll[3] += tz; // AT
    tz = n[4] * qp[4].ln_k[0];
    ll[4] += tz;  // CC
    ll[6] += tz;  // CT
    ll[9] += tz; // TT
  }
  if (n[5]) {
    ll[4] += log(Z[0] + qp[5].k) * n[5]; // CC
    double tz = log(0.5 * Z[2] + qp[5].k) * n[5];
    ll[1] += tz;                              // AC
    ll[5] += tz;                              // CG
    ll[6] += log(0.5 * Z[1] + qp[5].k) * n[5]; // CT
    tz = n[5] * qp[5].ln_k[0];
    ll[0] += tz;  // AA
    ll[2] += tz;  // AG
    ll[3] += tz;  // AT
    ll[7] += tz;  // GG
    ll[8] += tz;  // GT
    ll[9] += tz; // TT
  }
  if (n[6]) {
    ll[7] += log(Z[3] + qp[6].k) * n[6]; // GG
    double tz = log(0.5 * Z[5] + qp[6].k) * n[6];
    ll[5] += tz; // CG
    ll[8] += tz; // TG
    ll[2] += log(0.5 * Z[4] + qp[6].k) * n[6]; // AG
    tz = n[6] * qp[6].ln_k[0];
    ll[0] += tz;  // AA
    ll[1] += tz;  // AC
    ll[3] += tz;  // AT
    ll[4] += tz;  // CC
    ll[6] += tz;  // CT
    ll[9] += tz; // TT
  }
  if (n[7]) {
    ll[9] += n[7] * qp[7].ln_k[2];  // TT
    ll[4] += log(1.0 - Z[0] + qp[7].k) * n[7]; // CC
    ll[6] += log(1.0 - 0.5 * Z[1] + qp[7].k) * n[7]; // CT
    double tz = log(0.5 * (1.0 - Z[2]) + qp[7].k) * n[7];
    ll[1] += tz; // AC
    ll[5] += tz; // CG
    tz = n[7] * qp[7].ln_k[1];
    ll[3] += tz; // AT
    ll[8] += tz; // GT
    tz = n[7] * qp[7].ln_k[0];
    ll[0] += tz; // AA
    ll[2] += tz; // AG
    ll[7] += tz; // GG
  }
  double max = ll[0];
  int mx = 0;
  for (int i = 1; i < 10; i++) {
    if (ll[i] > max) {
      max = ll[i];
      mx = i;
    }
  }
  gt->max_gt = mx;
  double sum = 0.0;
  for (int i = 0; i < 10; i++) {
    sum += exp(ll[i] - max);
  }
  sum = log(sum);
  gt->sum = sum + max;
  for (int i = 0; i < 10; i++) {
    gt->gt_prob[i] = (ll[i] - max - sum);
  }
}

static int gt_idx[10][2] = {
  {-1, -1}, // AA
  {2, -1}, // AC
  {-1, 1}, // AG
  {-1, -1}, // AT
  {0, -1}, // CC
  {2, 2}, // CG
  {1, -1}, // CT
  {-1, 0}, // GG
  {-1, 2}, // GT
  {-1, -1} // TT
};

double get_meth(gt_meth *g, int idx) {
  double m = -1.0;
  int i = gt_idx[g->max_gt][idx];
  if(i >= 0) m = idx ? g->gmeth[i] : g->cmeth[i];
  return m;
}

// Calculate combined methylation for a CpG using information from both strands
// if available, taking account of the called genotypes.  If information is not
// available from both strands, use the single site estimate of methylation
void calc_cpg_meth(args_t *args, int ns, cpg_prob *cpg, gt_meth *g1, gt_meth *g2) {
  double wval[3] = {1.0, 1.0, 0.5};
  double pval[3] = {1.0, 0.5, 1.0};
  for(int ix = 0; ix < ns; ix++) {
    if(g1[ix].skip || g2[ix].skip) continue;
    int gt1 = g1[ix].max_gt;
    int gt2 = g2[ix].max_gt;
    cpg[ix].max_gt[0] = gt1;
    cpg[ix].max_gt[1] = gt2;
    cpg[ix].prob_best = g1[ix].gt_prob[gt1] + g2[ix].gt_prob[gt2];
    cpg[ix].prob_cg = g1[ix].gt_prob[4] + g2[ix].gt_prob[7];
    // Calc meth
    double n1[8], n2[8];
    qual_prob qp1[8], qp2[8];
    for (int i = 0; i < 8; i++) {
      n1[i] = (double)g1[ix].counts[i];
      n2[i] = (double)g2[ix].counts[i];
      qp1[i] = q_prob[g1[ix].aqual[i]];
      qp2[i] = q_prob[g2[ix].aqual[i]];
    }
    double l = 1.0 - args->under_conv;
    double t = args->over_conv;
    double g = (l - t) * 0.5;
    double f = (2.0 - l - t) * 0.5;
    double kc = qp1[5].k;
    double kt = qp1[7].k;
    double kg = qp2[6].k;
    double ka = qp2[4].k;
    int ix1 = gt_idx[gt1][0];
    int ix2 = gt_idx[gt2][1];
    if(ix1 >= 0) {
      double w1 = wval[ix1];
      double p = pval[ix1];
      if(ix2 >= 0) {
	double w2 = wval[ix2];
	double q = pval[ix2];
	// Get initial estimate
	double m1 = (n1[5] + n2[6]) / (n1[5] + n2[6] + n1[7] * p + n2[4] * q);
	double m = asin(2.0 * m1 - 1.0);
	// Maximize using NR
	for(int it = 0; it < 100; it++) {
	  double cosm = cos(m);
	  double sinm = sin(m);
	  double A = f + g * sinm;
	  double nm1 = g * p * w1 * cosm;
	  double d1 = p * w1 * A + kc;
	  double d2 = w1 * (1.0 - p * A) + kt;
	  double nm3  = g * q * w2 * cosm;
	  double d3 = q * w2 * A + kg;
	  double d4 = w2 * (1.0 - q * A) + ka;
	  double grad = nm1 * (n1[5] / d1 - n1[7] / d2) + nm3 * (n2[6] / d3 - n2[4] / d4);
	  if(fabs(grad) < 1.0e-8) {
	    m1 = 0.5 * (sinm + 1.0);
	    break;
	  }
	  double h = n1[5] * (nm1 * nm1 / d1 + g * p * w1 * sinm) / d1 + n1[7] * (nm1 * nm1 / d2 - g * p * w1 * sinm) / d2 +
	    n2[6] * (nm3 * nm3 / d3 + g * q * w2 * sinm) / d3 + n2[4] * (nm3 * nm3 / d4 - g * q * w2 * sinm) / d4;
	  m += grad / h;
	}
	cpg[ix].m = m1;
      } else {
	// Only the C+ has an estimate of methylation
	double m1 = g1->cmeth[ix1];
	cpg[ix].m = m1;
      }
    } else if(ix2 >= 0) {
	// Only the C- has an estimate of methylation
	double m1 = g2->cmeth[ix2];
	cpg[ix].m = m1;
    } else {
      // No valud esetimates on either strand
      cpg[ix].m = -1.0;
    }
  }
}

