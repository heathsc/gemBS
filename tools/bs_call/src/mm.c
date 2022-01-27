/*
 * mm.c
 *
 *  Created on: Feb 7, 2020
 *      Author: heath
 */

#include <stdlib.h>
#include <stdio.h>

void *safe_malloc(size_t sz, int line, char *file, const char *func) {
	void *p = malloc(sz);
	if(!p) {
		fprintf(stderr, "malloc() failed.  Requested %zu bytes from %s (%s:%d)\n", sz, func, file, line);
		exit(-1);
	}
	return p;
}

void *safe_calloc(size_t nm, size_t sz, int line, char *file, const char *func) {
	void *p = calloc(nm, sz);
	if(!p) {
		fprintf(stderr, "calloc() failed.  Requested %zu bytes from %s (%s:%d)\n", nm * sz, func, file, line);
		exit(-1);
	}
	return p;
}

void *safe_realloc(void *ptr, size_t sz, int line, char *file, const char *func) {
	void *p = realloc(ptr, sz);
	if(!p) {
		fprintf(stderr, "realloc() failed.  Requested %zu byte from %s (%s:%d)\n", sz, func, file, line);
		free(ptr);
		exit(-1);
	}
	return p;
}
