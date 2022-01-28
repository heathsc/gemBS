/*
 * stats_utils.c
 *
 *  Created on: Nov 24, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <math.h>

#include "gem_tools.h"
#include "bs_call.h"

void lfact_store_init(double * const lfact_store) {
	lfact_store[0] = lfact_store[1] = 0.0;
	double l = 0.0;
	for (int i = 2; i < LFACT_STORE_SIZE; i++) {
		l += log((double)i);
		lfact_store[i] = l;
	}
}

#define lfact(x) lfact2(x,lfact_store)

double fisher(int *const c, const double * const lfact_store) {
	int row[2], col[2];
	row[0] = c[0] + c[1];
	row[1] = c[2] + c[3];
	col[0] = c[0] + c[2];
	col[1] = c[1] + c[3];
	int n = row[0] + row[1];
	if(n == 0) return 1.0;
	double delta = (double)c[0] - (double)(row[0] * col[0]) / (double)n;
	double knst = lfact(col[0]) + lfact(col[1]) + lfact(row[0]) + lfact(row[1]) - lfact(n);
	double l = exp(knst - lfact(c[0]) - lfact(c[1]) - lfact(c[2]) - lfact(c[3]));
	double p = l;
	if (delta > 0.0) {
		// Decrease counter diagonal elements until zero (this will increase delta)
		int mn = c[1] < c[2] ? c[1] : c[2];
		for (int i = 0; i < mn; i++) {
			l *= (double)((c[1] - i) * (c[2] - i)) /
					(double)((c[0] + i + 1) * (c[3] + i + 1));
			p += l;
		}
		mn = c[0] < c[3] ? c[0] : c[3];
		// Calculate amount required to increase delta by decreasing leading
		// diagonal elements
		int k = ceil(2.0 * delta);
		if (k <= mn) {
			c[0] -= k;
			c[3] -= k;
			c[1] += k;
			c[2] += k;
			l = exp(knst - lfact(c[0]) - lfact(c[1]) - lfact(c[2]) - lfact(c[3]));
			p += l;
			for (int i = 0; i < mn - k; i++) {
				l *= (double)((c[0] - i) * (c[3] - i)) /
						(double)((c[1] + i + 1) * (c[2] + i + 1));
				p += l;
			}
		}
	} else {
		// Decrease leading diagonal elements until zero (this will increase delta)
		int mn = c[0] < c[3] ? c[0] : c[3];
		for (int i = 0; i < mn; i++) {
			l *= (double)((c[0] - i) * (c[3] - i)) /
					(double)((c[1] + i + 1) * (c[2] + i + 1));
			p += l;
		}
		mn = c[1] < c[2] ? c[1] : c[2];
		// Calculate amount required to increase delta by decreasing counter
		// diagonal elements
		int k = ceil(-2.0 * delta);
		if (!k)
			k = 1;
		if (k <= mn) {
			c[0] += k;
			c[3] += k;
			c[1] -= k;
			c[2] -= k;
			l = exp(knst - lfact(c[0]) - lfact(c[1]) - lfact(c[2]) - lfact(c[3]));
			p += l;
			for (int i = 0; i < mn - k; i++) {
				l *= (double)((c[1] - i) * (c[2] - i)) /
						(double)((c[0] + i + 1) * (c[3] + i + 1));
				p += l;
			}
		}
	}
	return p;
}

