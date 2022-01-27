/*
 * dbSNP_buffer.c
 *
 *  Created on: Feb 7, 2020
 *      Author: heath
 */

#include "dbSNP_idx.h"
#include "dbSNP_utils.h"

#ifndef kroundup_size_t
#define kroundup_size_t(x) (--(x),                    \
		(x)|=(x)>>(sizeof(size_t)/8), /*  0 or  1 */ \
		(x)|=(x)>>(sizeof(size_t)/4), /*  1 or  2 */ \
		(x)|=(x)>>(sizeof(size_t)/2), /*  2 or  4 */ \
		(x)|=(x)>>(sizeof(size_t)),   /*  4 or  8 */ \
		(x)|=(x)>>(sizeof(size_t)*2), /*  8 or 16 */ \
		(x)|=(x)>>(sizeof(size_t)*4), /* 16 or 32 */ \
		++(x))
#endif

void init_buffer(buffer_t * const buf, size_t sz) {
	kroundup_size_t(sz);
	buf->size = sz;
	buf->mem = malloc(sz);
	buf->len = 0;
}

buffer_t *new_buffer(size_t sz) {
	if(sz < 32) sz = 32;
	buffer_t *buf = malloc(sizeof(buffer_t));
	init_buffer(buf, sz);
	return buf;
}

void destroy_buffer(buffer_t * const buf) {
	free(buf->mem);
	free(buf);
}

void resize_buffer(buffer_t * const buf, size_t sz) {
	if(sz > buf->size) {
		kroundup_size_t(sz);
		buf->mem = realloc(buf->mem, sz);
		buf->size = sz;
	}
}

void swap_buffers(buffer_t * const buf1, buffer_t * const buf2) {
	buffer_t t;
	memcpy(&t, buf1, sizeof(buffer_t));
	memcpy(buf1, buf2, sizeof(buffer_t));
	memcpy(buf2, &t, sizeof(buffer_t));
}
