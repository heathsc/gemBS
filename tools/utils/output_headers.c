/*
 * output_headers.c
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

static char *copy_and_strip_quotes(char *s) {
	if(!s) return s;
	size_t l = strlen(s);
	if(l > 1) {
		if((s[0] == '\"' && s[l-1] =='\"') || (s[0] == '\'' && s[l-1] =='\'')) {
			(s++)[--l] = 0;
		}
	}
	char *s1 = malloc(l + 1);
	if(s1 != NULL) memcpy(s1, s, l + 1);
	return s1;
}

static void print_file_header(htsFile *fp, kstring_t *s, int ns, char **names) {
	if(fp != NULL) {
		kputs("Contig\tPos0\tPos1\tRef", ks_clear(s));
		for(int i = 0; i < ns; i++) {
			char *name = names[i];
			ksprintf(s, "\t%s:Call\t%s:Flags\t%s:Meth\t%s:non_conv\t%s:conv\t%s:support_call\t%s:total", name, name, name, name, name, name, name);
		}
		kputc('\n', s);
		ks_output(fp, s);
	}
}

static void print_bedmethyl_headers(args_t *args) {
	kstring_t *s = ks_clear(&args->pr_str[2]);
	if(args->bedmethyl_track_line == NULL) {
		char *sample_name = NULL;
		char *sample_desc = NULL;
		char *sample_bc = NULL;
		// Try and get sample info from VCF file headers
		bcf_hdr_t *h = args->hdr;
		for(int i = 0; i < h->nhrec; i++) {
			bcf_hrec_t *hr = h->hrec[i];
			if(hr->type == BCF_HL_STR) {
				if(!strcmp(hr->key, "bs_call_sample_info")) {
					int ix = bcf_hrec_find_key(hr, "ID");
					if(ix >= 0) {
						sample_bc = copy_and_strip_quotes(hr->vals[ix]);
						ix = bcf_hrec_find_key(hr, "SM");
						if(ix >= 0) sample_name = copy_and_strip_quotes(hr->vals[ix]);
						ix = bcf_hrec_find_key(hr, "DS");
						if(ix >= 0) sample_desc = copy_and_strip_quotes(hr->vals[ix]);
					}
				}
			}
		}
		if(sample_name == NULL) sample_name = strdup(h->samples[0]);
		if(sample_desc == NULL) sample_desc = strdup(sample_name);
		ksprintf(s, "track name=\"%s\" description=\"%s\" visibility=2 itemRgb=\"On\"\n", sample_desc, sample_name);
		if(sample_bc) free(sample_bc);
		if(sample_name) free(sample_name);
		args->bedmethyl_desc = sample_desc;
	} else {
		char *line = args->bedmethyl_track_line;
		size_t l = strlen(line);
		if(l > 1 && line[l - 1] == '\n') line[--l] = 0;
		if(!strncmp(line, "track ", 6)) line += 6;
		ksprintf(s, "track %s\n", line);
	}
	for(bedmethyl_type t = BEDMETHYL_CPG; t <= BEDMETHYL_CHH; t++) {
		htsFile *fp = args->bedmethylfiles[t - 1];
		if(fp != NULL) ks_output(fp, s);
	}
}

void print_headers(args_t *args) {
	int ns = bcf_hdr_nsamples(args->hdr);
	if(args->cpgfile) print_file_header(args->cpgfile, &args->pr_str[0], ns, args->hdr->samples);
	if(args->noncpgfile) print_file_header(args->noncpgfile, &args->pr_str[1], ns, args->hdr->samples);
	if(args->bedmethyl) print_bedmethyl_headers(args);
}

