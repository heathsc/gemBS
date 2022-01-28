#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <sys/wait.h>

#include "utils.h"

static void cat_file(FILE *fp) {
	while(!feof(fp)) {
		int c = fgetc(fp);
		if(c == EOF) break;
		putchar(c);
	}
}

int main(int argc, char *argv[]) {
	int err = 0;

	for(int ix = 1; ix <= argc; ix++) {
		if(ix == argc) {
			if(argc == 1) {
				cat_file(stdin);
			}
		} else {
			bool flag;
			FILE *fp = open_readfile(argv[ix], &flag);
			if(fp == NULL) {
				err = errno;
				break;
			}
			cat_file(fp);
			fclose(fp);
			if(flag) {
				while(waitpid(-1, NULL, 0) > 0);
			}
		}
	}
	return err;
}
