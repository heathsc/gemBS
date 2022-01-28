/*
 * files.c
 *
 *  Created on: Jan 24, 2020
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

#include "snpxtr.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"

htsFile *open_ofile(char ** const name, bool compress, sargs_t * const a) {

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



