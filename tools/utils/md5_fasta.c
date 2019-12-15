/*
 * md5_fasta.c
 *
 *  Created on: Dec 12, 2019
 *      Author: heath
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <errno.h>
#include <getopt.h>
#include <openssl/md5.h>

#include "utils.h"

#ifndef PATH_MAX
#define PATH_MAX FILENAME_MAX
#endif

static void usage(FILE * const f) {
	fputs("usage:\n md5_fasta [option]... [file]...\n", f);
	fputs("  -o|--output          <output file>\n", f);
	fputs("  -s|--stream          stream input to stdout\n", f);
	fputs("  -p|--populate-cache  populate local cache\n", f);
	fputs("  -h|--help            print this message\n", f);
}

static char btab[256];
static char path[PATH_MAX];

static int get_cache_fname(char *path, const char * cache, const char * str) {
	size_t s = PATH_MAX;
	size_t ln;
	char *init = path;
	char *p;
	while((p = strchr(cache, '%'))) {
		ln = p - cache;
		if(ln >= s) return -1;
		while(cache < p) *path++ = *cache++;
		s -= ln;
		if(*++p == 's') {
			ln = strlen(str);
			if(ln > s) return -1;
			strcpy(path, str);
			s -= ln;
			path += ln;
			str += ln;
			p++;
		} else if(*p >= '0' && *p <= '9') {
			char *p1;
			long l = strtol(p, &p1, 10);
			ln = strlen(str);
			if(l > ln) l = ln;
			if(*p1 == 's') {
				if(l >= s) return -1;
				s -= l;
				while(l--) *path++ = *str++;
				p = p1 + 1;
			} else {
				if(s < 3) return -1;
				*path++ = '%';
				*path++ = *p++;
			}
		} else {
			if(s < 3) return -1;
			*path++ = '%';
			*path++ = *p++;
		}
		cache = p;
	}
	ln = strlen(cache);
	if(ln >= s) return -1;
	s -= ln;
	while(ln--) *path++ = *cache++;
	bool need_slash = (*str && path > init && path[-1] != '/');
	ln = strlen(str) + need_slash ? 1 : 0;
	if(ln >= s) return -1;
	if(need_slash) *path++ = '/';
	strcpy(path, str);
	return 0;
}

static int mk_path(char *newpath) {
	static char path[PATH_MAX];
	if(strlen(newpath) >= PATH_MAX) return -1;
	struct stat sb;
	char *p, *p1 = path;
	if(*newpath == '/') *p1++ = *newpath++;
    while((p = strchr(newpath, '/'))) {
    	while(newpath < p) *p1++ = *newpath++;
    	*p1 = 0;
    	if(stat(path, &sb) == 0) {
    		if(!S_ISDIR(sb.st_mode)) return -1;
    	} else {
    		if(errno != ENOENT) return -1;
    		if(mkdir(path, 0777)) return -1;
    	}
    	*p1++ = *newpath++;
    }
	return 0;
}

static void handle_md5(const char * const ctg, MD5_CTX * const ctx, size_t tlen,
		const char * const cache, const char * const ref_buf, size_t ref_len, FILE *fout) {
	unsigned char b[16];
	char md5[33];
	char *hex_digits="0123456789abcdef";
	MD5_Final(b, ctx);
	const char *p = ctg;
	int k = 0;
	for(int i = 0; i < 16; i++) {
		md5[k++] = hex_digits[b[i] >> 4];
		md5[k++] = hex_digits[b[i] & 0xf];
	}
	md5[k] = 0;
	while(*p && !isspace(*p)) fputc(*p++, fout);
	fprintf(fout,"\tLN:%zu\tM5:%s", tlen, md5);
	while(*p) {
		while(*p && isspace(*p)) p++;
		if(!(strncmp(p, "AS:", 3) && strncmp(p, "SP:", 3))) {
			fputc('\t', fout);
			while(*p && !isspace(*p)) fputc(*p++, fout);
		} else while(*p && !isspace(*p)) p++;
	}
	fputc('\n', fout);
	if(cache && ref_buf) {
		if(get_cache_fname(path, cache, md5) == 0) {
			struct stat sb;
			if(stat(path, &sb) != 0) {
				if(errno == ENOENT) {
					fprintf(stderr, "Creating cache file %s\n", path);
					if(mk_path(path) == 0) {
						FILE *fp = fopen(path, "w");
						if(fp) {
							size_t s = fwrite(ref_buf, 1, ref_len, fp);
							fclose(fp);
							if(s != ref_len) {
								fprintf(stderr, "Error writing cache file %s: %s\n", path, strerror(errno));
								unlink(path);
							}
						} else fprintf(stderr, "Could not create cache file %s: %s\n", path, strerror(errno));
					} else {
						fprintf(stderr, "Could not create directory path for cache file %s\n", path);
					}
				} else {
					fprintf(stderr, "Could not access directory path for cache file %s: %s\n", path, strerror(errno));
				}
			}
		} else {
			fprintf(stderr, "Could not create directory path for cache file %s: path too long\n", path);
		}
	}
}

static void process_file(FILE * const fp, FILE * const fout, const bool stream, const char * const cache) {
	char *buf = NULL;
	char *ref_buf = NULL;
	size_t ref_buf_size = 0, ref_len = 0;
	if(cache) {
		ref_buf_size = 16384;
		ref_buf = malloc(ref_buf_size);
		if(!ref_buf) {
			fprintf(stderr, "md5_fasta:process_file() Out of memory");
			exit(-1);
		}
	}
	size_t buf_size = 0, tlen = 0;
	ssize_t l;
	char *ctg = NULL;
	MD5_CTX ctx;
	// Process header lines - no conversion
	while(1) {
		l = getline(&buf, &buf_size, fp);
		if(l < 0) break;
		if(stream) fputs(buf, stdout);
		if(buf[0] == '>') {
			if(ctg) {
				handle_md5(ctg, &ctx, tlen, cache, ref_buf, ref_len, fout);
				ref_len = 0;
				free(ctg);
			}
			if(!buf[0]) continue;
			ctg = strdup(buf + 1);
			MD5_Init(&ctx);
			tlen = 0;
		} else {
			if(!ctg) {
				fprintf(stderr,"md5_fasta:no header found");
				exit(-1);
			}
			// First, strip characters not between 33 and 126, and convert to upper case.
			char *p = buf, *p1 = buf;
			while(*p) {
				char c = btab[(int)*p++];
				if(c) *p1++ = c;
			}
			size_t len = p1 - buf;
			if(len) {
				MD5_Update(&ctx, buf, len);
				tlen += len;
				if(ref_buf) {
					if(len + ref_len > ref_buf_size) {
						ref_buf_size = (ref_buf_size + len) * 1.5;
						ref_buf = realloc(ref_buf, ref_buf_size);
						if(!ref_buf) {
							fprintf(stderr, "md5_fasta:process_file() Out of memory");
							exit(-1);
						}
					}
					memcpy(ref_buf + ref_len, buf, len);
					ref_len += len;
				}
			}
		}
	}
	if(ctg) {
		handle_md5(ctg, &ctx, tlen, cache, ref_buf, ref_len, fout);
		free(ctg);
	}
	if(buf != NULL) free(buf);
	if(ref_buf != NULL) free(ref_buf);
}

static char *get_cache_path(void) {
	char *ref_cache = getenv("REF_CACHE");
	if(!ref_cache || *ref_cache == 0) {
		// This is the same logic used in htslib/cram/cram_io.c
		char *ext = NULL;
		char *base = getenv("XDG_CACHE_HOME");
		if(!(base && *base)) {
			base = getenv("HOME");
			if(base && *base) ext = "/.cache";
		}
		if(!(base && *base)) base = getenv("TMPDIR");
		if(!(base && *base)) base = getenv("TEMP");
		if(!(base && *base)) base = "/tmp";
		if(!ext) ext = "";
		ref_cache = malloc(PATH_MAX);
		snprintf(ref_cache, PATH_MAX, "%s%s/hts-ref/%%2s/%%2s/%%s", base, ext);
	}
	return ref_cache;
}

int main(int argc, char* argv[]) {
	struct option longopts[] = {
			{"output", required_argument, 0, 'o'},
			{"stream", no_argument, 0, 's'},
			{"populate-cache", no_argument, 0, 'p'},
			{"help", no_argument, 0, 'h'},
			{"usage", no_argument, 0, 'h'},
			{NULL, 0, 0, 0}
	};
	int err = 0;
	char *output = NULL;
	bool stream = false;
	bool populate_cache = false;
	for(int i = 33; i < 127; i++) btab[i] = toupper(i);
	int c;
	while(!err && (c = getopt_long(argc, argv, "o:psh?", longopts, 0)) != -1) {
		switch(c) {
		case 'o':
			output = optarg;
			break;
		case 's':
			stream = true;
			break;
		case 'p':
			populate_cache = true;
			break;
		case 'h':
		case '?':
			usage(stdout);
		}
	}
	char *cache = NULL;
	if(populate_cache) cache = get_cache_path();
	FILE *fout = NULL;
	if(output) fout = fopen(output, "w");
	if(fout == NULL) fout = stdout;
	for(int ix = optind; ix <= argc; ix++) {
		if(ix == argc) {
			if(argc == optind) process_file(stdin, fout, stream, cache);
		} else {
			bool flag;
			FILE *fp = open_readfile(argv[ix], &flag);
			if(fp == NULL) {
				err = errno;
				break;
			}
			process_file(fp, fout, stream, cache);
			fclose(fp);
			if(flag) while(waitpid(-1, NULL, 0) > 0);
		}
	}
	if(fout != stdout) fclose(fout);
	return err;
}
