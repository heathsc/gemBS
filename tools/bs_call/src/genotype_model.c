#include <stdio.h>
#include <math.h>
#include <getopt.h>

#include "gem_tools.h"
#include "bs_call.h"

static qual_prob q_prob[MAX_QUAL + 1];

void fill_base_prob_table(void) {
  for (int q = 0; q <= MAX_QUAL; q++) {
    double e = exp(-.1 * (double)q * LOG10);
    if(e > .5) e = .5;
    double k = e / (3.0 - 4.0 * e);
    q_prob[q].e = e;
    q_prob[q].k = k;
    q_prob[q].ln_k = log(k);
    q_prob[q].ln_k_half = log(0.5 + k);
    q_prob[q].ln_k_one = log(1.0 + k);
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

void calc_gt_prob(gt_meth *gt, const sr_param * const param, char rf) {
  qual_prob qp[8];
  for(int i = 0; i < 8; i++) qp[i] = q_prob[gt->qual[i]];
  double l = 1.0 - param->under_conv;
  double t = param->over_conv;

  /***********************************************************************************
   * Base and methylation frequencies are described by 5 parameters: w, p, q, mc, mg
   * 
   * Let n(X) be the count for base X, and N the total number of bases seen
   * w = (n(C) + n(T)) / N
   * p = n(C) / (n(C) + n(T))
   * q = n(G) / (n(A) + n(G))
   * mc is the proportion of methylated Cs on the top strand
   * mg is the proportion of methylated Cs on the bottom strand
   *
   * Base frequencies are therefore:
   *  f(A) = (1 - w) * (1 - q)
   *  f(C) = w * p
   *  f(G) = (1 - w) * q
   *  f(T) = w * (1 - p)
   *
   * All 5 parameters are ratios are are therefore independently constrained 
   * to be between 0 and 1.
   * 
   * We first maximize the full model, allowing w, p, q, mc and mg to take 
   * any legal value.  The log likelihood of this model is l_full.
   *
   * We then calculate the marginal likelihood for the 10 possible genetic models compatible
   * with a diploid state (thereby fixing w, p, q) and maximizing the likelihood over (mc, mg).
   *
   * The called genotype is that with the highest likelihood 
   * The phred score is calculated as the phred scaled posterior genotype probability (considering
   * only the 10 possible diploid genotypes)
   *
   ***********************************************************************************/

  double ref_bias = param->ref_bias;
  double n[8];
  for (int i = 0; i < 8; i++) n[i] = (double)gt->counts[i];
  double ll[10];
  memset(ll, 0, sizeof(double) * 10);
  // Add in prior from reference
  {
    double lrb = log(ref_bias);
    double lrb1 = log(0.5 * (1.0 + ref_bias));
    switch (rf) {
    case 1: // A
      ll[0] = lrb;
      ll[1] = ll[2] = ll[3] = lrb1;
      break;
    case 2: // C
      ll[4] = lrb;
      ll[1] = ll[5] = ll[6] = lrb1;
      break;
    case 3: // G
      ll[7] = lrb;
      ll[2] = ll[5] = ll[8] = lrb1;
      break;
    case 4: // T
      ll[9] = lrb;
      ll[3] = ll[6] = ll[8] = lrb1;
      break;
    }
  }
  if (n[0]) {
    ll[0] += n[0] * qp[0].ln_k_one; // AA
    double tz = n[0] * qp[0].ln_k_half; 
    ll[1] += tz; // AC
    ll[2] += tz; // AG
    ll[3] += tz; // AT
    tz = n[0] * qp[0].ln_k;
    ll[4] += tz;  // CC
    ll[5] += tz;  // CG
    ll[6] += tz;  // CT
    ll[7] += tz;  // GG
    ll[8] += tz;  // GT
    ll[9] += tz; // TT
  }
  if (n[1]) {
    ll[4] += n[1] * qp[1].ln_k_one; // CC
    double tz = n[1] * qp[1].ln_k_half;
    ll[1] += tz; // AC
    ll[5] += tz; // CG
    ll[6] += tz; // CT
    tz = n[1] * qp[1].ln_k;
    ll[0] += tz;  // AA
    ll[2] += tz;  // AG
    ll[3] += tz;  // AT
    ll[7] += tz;  // GG
    ll[8] += tz;  // GT
    ll[9] += tz; // TT
  }
  if (n[2]) {
    ll[7] += n[2] * qp[2].ln_k_one; // GG
    double tz = n[2] * qp[2].ln_k_half;
    ll[2] += tz; // AG
    ll[5] += tz; // CG
    ll[8] += tz; // TG
    tz = n[2] * qp[2].ln_k;
    ll[0] += tz;  // AA
    ll[1] += tz;  // AC
    ll[3] += tz;  // AT
    ll[4] += tz;  // CC
    ll[6] += tz;  // CT
    ll[9] += tz; // TT
  }
  if (n[3]) {
    ll[9] += n[3] * qp[3].ln_k_one; // TT
    double tz = n[3] * qp[3].ln_k_half;
    ll[3] += tz; // AT
    ll[6] += tz; // CT
    ll[8] += tz; // GT
    tz = n[3] * qp[3].ln_k;
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
  }
  if (n[4] + n[6] > 0.0) {
    get_Z(n[6], n[4], qp[6].k, qp[4].k, l, t, Z+3);
  }

  if (n[4]) {
    ll[0] += n[4] * qp[4].ln_k_one; // AA
    ll[2] += log(1.0 - 0.5 * Z[4] + qp[4].k) * n[4]; // AG
    ll[7] += log(1.0 - Z[3] + qp[4].k) * n[4]; // GG
    double tz = log(0.5 * (1.0 - Z[5]) + qp[4].k) * n[4];
    ll[5] += tz; // CG
    ll[8] += tz; // GT
    tz = n[4] * qp[4].ln_k_half;
    ll[1] += tz; // AC
    ll[3] += tz; // AT
    tz = n[4] * qp[4].ln_k;
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
    tz = n[5] * qp[5].ln_k;
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
    tz = n[6] * qp[6].ln_k;
    ll[0] += tz;  // AA
    ll[1] += tz;  // AC
    ll[3] += tz;  // AT
    ll[4] += tz;  // CC
    ll[6] += tz;  // CT
    ll[9] += tz; // TT
  }
  if (n[7]) {
    ll[9] += n[7] * qp[7].ln_k_one;  // TT
    ll[4] += log(1.0 - Z[0] + qp[7].k) * n[7]; // CC
    ll[6] += log(1.0 - 0.5 * Z[1] + qp[7].k) * n[7]; // CT
    double tz = log(0.5 * (1.0 - Z[2]) + qp[7].k) * n[7];
    ll[1] += tz; // AC
    ll[5] += tz; // CG
    tz = n[7] * qp[7].ln_k_half;
    ll[3] += tz; // AT
    ll[8] += tz; // GT
    tz = n[7] * qp[7].ln_k;
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
  for (int i = 0; i < 10; i++) sum += exp(ll[i] - max);
  sum = log(sum);
  for (int i = 0; i < 10; i++) {
    gt->gt_prob[i] = (ll[i] - max - sum) / LOG10;
  }
}


