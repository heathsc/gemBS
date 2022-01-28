/*
 * stats.c
 *
 *  Created on: Dec 26, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "mextr.h"

void init_stats(args_t *a) {
	a->stats = calloc((size_t)1, sizeof(stats_t));
}

void write_stats(args_t *a) {
	if(a->stats != NULL) {
		a->reportfile = a->reportfilename == NULL ? NULL : fopen(a->reportfilename, "w");
		if(a->reportfile != NULL) {
			FILE *fp = a->reportfile;
			stats_t *st = a->stats;
			fprintf(fp,"{\n\t\"TotalSites\": %" PRIu64 ",\n", st->n_sites);
			fprintf(fp,"\t\"SitesPassed\": %" PRIu64 "\n", st->n_sites_pass);
			fputs("}\n", fp);
			fclose(fp);
		}
	}
}

