/*
 * output_utils.c
 *
 *  Created on: Dec 26, 2019
 *      Author: heath
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "mextr.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"

int calc_phred(double z) {
	int phred;
	if(z <= 0.0) phred = 255;
	else {
		phred = (int)(-10.0 * log(z) / LOG10);
		if(phred > 255) phred = 255;
	}
	return phred;
}

double *get_prob_dist(int ns, double *Q[]) {
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

char trans_base[256] = {
		['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A',
		['Y'] = 'R', ['R'] = 'Y', ['S'] = 'S', ['W'] = 'W', ['K'] = 'M', ['M'] = 'K',
		['B'] = 'V', ['V'] = 'B', ['D'] = 'H', ['H'] = 'D', ['N'] = 'N', ['.'] = '.'
};

