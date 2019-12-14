#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/wait.h>
#include <getopt.h>

#include "utils.h"
#include "uthash.h"

// Strip illegal characters from Read IDs in SAM file
// Valid characters are [!-?A-~]

// Option to edit SAM headers, adding extra information to the @SQ lines

#define NUM_SQTAGS 4
static char *sqtags[NUM_SQTAGS] = {
		"LN", "M5", "AS", "SP"
};

typedef struct {
	char *name;
	char *tags[NUM_SQTAGS];
	UT_hash_handle hh;
} ctg_t;

ctg_t *process_ctg_file(char *name) {
	ctg_t *ctgs = NULL;
	bool flag;
	FILE *fp = open_readfile(name, &flag);
	if(fp == NULL) {
		fprintf(stderr, "Could not open %s for reading\n", name);
		exit(-1);
	}
	char *buf = NULL;
	size_t buf_size = 0, tlen = 0;
	ssize_t l;
	tokens *tok = NULL;
	while(1) {
		l = getline(&buf, &buf_size, fp);
		if(l < 0) break;
		tok = tokenize(buf, '\t', tok);
		if(tok->n_tok > 1) {
			ctg_t *ct = NULL;
			HASH_FIND_STR(ctgs, tok->toks[0], ct);
			if(ct != NULL) {
				fprintf(stderr, "process_ctg_file(): error - duplicate contig %s\n", tok->toks[0]);
				exit(-1);
			}
			ct = malloc(sizeof(ctg_t));
			ct->name = strdup(tok->toks[0]);
			for(int i = 0; i < NUM_SQTAGS; i++) ct->tags[i] = NULL;
			for(int i = 1; i < tok->n_tok; i++) {
				const char * const p = tok->toks[i];
				for(int j = 0; j < NUM_SQTAGS; j++) {
					if(!strncmp(p, sqtags[j], 2) && p[2] == ':') {
						ct->tags[j] = strdup(p + 3);
						break;
					}
				}
			}
			HASH_ADD_KEYPTR(hh, ctgs, ct->name, strlen(ct->name), ct);
		}
	}
	fclose(fp);
	if(flag) while(waitpid(-1, NULL, 0) > 0);
	if(tok != NULL) free_tokens(tok);
	if(buf != NULL) free(buf);
	return ctgs;
}

int main(int argc, char *argv[]) {
	FILE *fp = stdin;
	char *buf = NULL;
	size_t buf_size = 0;
	ssize_t l;
	ctg_t *ctgs = NULL;

	if(argc > 1) ctgs = process_ctg_file(argv[1]);
	// Process header lines - no conversion
	while(1) {
		l = getline(&buf, &buf_size, fp);
		if(l < 0) return 0;
		if(buf[0] != '@') break;
		bool pflag = true;
		if(l > 8 && !strncmp(buf + 1, "SQ\tSN:", 6)) {
			char *p = buf + 7;
			char *p1 = p;
			while(*p1 && *p1 != '\t' && *p1 != '\n') p1++;
			size_t l = p1 - p;
			ctg_t *ct = NULL;
			HASH_FIND(hh, ctgs, p, l, ct);
			if(ct) {
				pflag = false;
				int mask = 0;
				char c = *p1;
				*p1 = 0;
				fputs(buf, stdout);
				while(c == '\t') {
					p1++;
					p = p1;
					for(int j = 0; j < NUM_SQTAGS; j++) {
						if(!strncmp(p1, sqtags[j], 2) && p1[2] == ':') {
							mask |= (1 << j);
							break;
						}
					}
					while(*p1 && *p1 != '\t' && *p1 != '\n') p1++;
					c = *p1;
					*p1 = 0;
					printf("\t%s", p);
				}
				int j = 0;
				for(int j = 0; j < NUM_SQTAGS; j++) {
					if(ct->tags[j] != NULL && !(mask & (1 << j))) printf("\t%s:%s", sqtags[j], ct->tags[j]);
				}
				fputc('\n', stdout);
			}
		}
		if(pflag) fputs(buf, stdout);
	}
	// Process the rest of the file
	while(l >= 0) {
		int i;
		bool found = false;
		for(i = 0; i < l && buf[i] != '\t'; i++) if((found = (buf[i] == '@' || buf[i] < '!' || buf[i] > '~'))) break;
		if(found) {
			int j = i;
			for(i = i + 1; i < l && buf[i] != '\t'; i++) {
				if(buf[i] != '@' && buf[i] >= '!' && buf[i] <= '~') buf[j++] = buf[i];
			}
			for(; i <= l; i++) buf[j++] = buf[i];
		}
		fputs(buf, stdout);
		l = getline(&buf, &buf_size, fp);
	}
	if(buf) free(buf);
	return 0;
}
