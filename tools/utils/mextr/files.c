/*
 * files.c
 *
 *  Created on: Dec 26, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <sys/wait.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <openssl/md5.h>

#include "mextr.h"
#include "bbi.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"

static htsFile *open_ofile(char ** const name, bool compress, args_t * const a) {

	htsFile *fp = NULL;
	if(name != NULL) {
		char *tname = *name;
		char mode[3] = {'w', 0, 0};
		bool stream = !strcmp(tname, "-");

		// If output is to a file then if compression has been asked for
		// we add '.gz' to the filename unless already present.  If compression
		// had not been asked for but the filename ends in '.gz' then we
		// turn on compression

		if(!stream) {
			// Check if file name ends in '.gz'
			char *p = strrchr(tname, '.');
			bool has_gz = p && !strcmp(p + 1, "gz");
			if(compress) {
				if(!has_gz) {
					tname = malloc(strlen(*name) + 4);
					sprintf(tname, "%s.gz", *name);
				}
			} else compress = has_gz;
		} else {
			// Turn off compression if output is to a terminal
			if(compress && isatty(fileno(stdout))) compress = false;
		}
		if(compress) mode[1] = 'z';
		hFILE *hfile = hopen(tname, mode);
		if(!hfile) error("Couldn't open output file: %s\n", stream ? "<stdout>" : tname);
		fp = hts_hopen(hfile, tname, mode);
		if(a->threads > 0) hts_set_opt(fp, HTS_OPT_THREAD_POOL, a->sr->p);
		if(tname != *name) *name = tname;
	}
	return fp;
}

void init_files(args_t * const a) {
	if(a->cpgfilename) a->cpgfile = open_ofile(&a->cpgfilename, a->compress, a);
	if(a->noncpgfilename) a->noncpgfile = open_ofile(&a->noncpgfilename, a->compress, a);
	if(a->bedmethyl != NULL) {
		if(!strcmp(a->bedmethyl, "-")) error("Bedmethyl files can not be output to stdout");
		bool compress = a->compress;
		char *p = strrchr(a->bedmethyl, '.');
		if(p && !strcmp(p + 1, "gz")) {
			*p = 0;
			compress = true;
		}
		p = strrchr(a->bedmethyl, '.');
		if(p && !strcmp(p + 1, "bed")) *p = 0;
		p = strrchr(a->bedmethyl, '_');
		if(p && !(strcmp(p + 1, "cpg") || strcmp(p + 1, "chg") || strcmp(p + 1, "chh")))*p = 0;
		const size_t l = strlen(a->bedmethyl) + 9 + (compress ? 3 : 0);
		const int ct = a->strand_specific ? 8 : 7;
		p = malloc(l * ct);
		a->bedmethylnames[0] = p;
		a->bedmethylnames[1] = p + l;
		a->bedmethylnames[2] = p + 2 * l;
		a->bigbednames[0] = p + 3 * l;
		a->bigbednames[1] = p + 4 * l;
		a->bigbednames[2] = p + 5 * l;
		a->bigwignames[0] = p + 6 * l;
		if(a->strand_specific) a->bigwignames[1] = p + 7 * l;

		sprintf(a->bedmethylnames[BEDMETHYL_CPG - 1], "%s_cpg.%s", a->bedmethyl, compress ? "bed.gz" : "bed");
		sprintf(a->bedmethylnames[BEDMETHYL_CHG - 1], "%s_chg.%s", a->bedmethyl, compress ? "bed.gz" : "bed");
		sprintf(a->bedmethylnames[BEDMETHYL_CHH - 1], "%s_chh.%s", a->bedmethyl, compress ? "bed.gz" : "bed");
		sprintf(a->bigbednames[BEDMETHYL_CPG - 1], "%s_cpg.bb", a->bedmethyl);
		sprintf(a->bigbednames[BEDMETHYL_CHG - 1], "%s_chg.bb", a->bedmethyl);
		sprintf(a->bigbednames[BEDMETHYL_CHH - 1], "%s_chh.bb", a->bedmethyl);
		if(a->strand_specific) {
			sprintf(a->bigwignames[0], "%s_pos.bw", a->bedmethyl);
			sprintf(a->bigwignames[1], "%s_neg.bw", a->bedmethyl);
		} else sprintf(a->bigwignames[0], "%s.bw", a->bedmethyl);
		for(int i = 0; i < 3; i++) {
			a->bedmethylfiles[i] = open_ofile(&a->bedmethylnames[i], compress, a);
			a->bigbedfiles[i]= fopen(a->bigbednames[i], "wb");
		}
		for(int i = 0; i < 2; i++) {
			if(a->bigwignames[i]) a->bigwigfiles[i] = fopen(a->bigwignames[i], "w");
		}
		init_bbi_header(a, true); // bigBed headers
		init_bbi_header(a, false); // bigWig headers
		init_cblocks(a, 4 * a->compress_threads);
	}
}

void close_files(args_t *a) {
	if(a->cpgfile != NULL) hts_close(a->cpgfile);
	if(a->noncpgfile != NULL) hts_close(a->noncpgfile);
	for(int i = 0; i < 3; i++) {
		if(a->bedmethylfiles[i] != NULL) hts_close(a->bedmethylfiles[i]);
		if(a->bigbedfiles[i] != NULL) fclose(a->bigbedfiles[i]);
	}
	for(int i = 0; i < 2; i++) {
		if(a->bigwigfiles[i] != NULL) fclose(a->bigwigfiles[i]);
	}
	while(waitpid(-1, NULL, 0) > 0);
}

#define MD5_BUF_SIZE 4096

void calc_stream_md5(FILE * const fp, char * const md5) {
	MD5_CTX ctx;
	MD5_Init(&ctx);
	uint8_t buf[MD5_BUF_SIZE];
	while(!feof(fp)) {
		size_t len = fread(buf, 1, MD5_BUF_SIZE, fp);
		if(len > 0) MD5_Update(&ctx, buf, len);
	}
	unsigned char b[16];
	const char *hex_digits="0123456789abcdef";
	MD5_Final(b, &ctx);
	int k = 0;
	for(int i = 0; i < 16; i++) {
		md5[k++] = hex_digits[b[i] >> 4];
		md5[k++] = hex_digits[b[i] & 0xf];
	}
	md5[k] = 0;
}

void calc_file_md5(char * const name) {
	 char *tname = malloc(strlen(name) + 5);
	 sprintf(tname, "%s.md5", name);
	 FILE *in = fopen(name, "rb");
	 FILE *out = NULL;
	 int err = 0;
	 if(in == NULL) {
		 fprintf(stderr, "calc_file_md5(): Could not open file %s for input: %s\n", name, strerror(errno));
		 err = 1;
	 } else {
		 out = fopen(tname, "wb");
		 if(out == NULL) {
			 fprintf(stderr, "calc_file_md5(): Could not open file %s for output: %s\n", tname, strerror(errno));
			 err = 2;
		 }
	 }
	 char md5[33];
	 if(!err) {
#ifndef __MACH__
	posix_fadvise(fileno(in), 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
		 calc_stream_md5(in, md5);
		 fprintf(out,"%s  %s\n", md5, name);
	 }
	 if(out) fclose(out);
	 if(in) fclose(in);
	 free(tname);
}

