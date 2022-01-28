/*
 * mm.h
 *
 *  Created on: Feb 7, 2020
 *      Author: heath
 */

#ifndef INCLUDE_MM_H_
#define INCLUDE_MM_H_

void *safe_malloc(size_t sz, int line, char *file, const char *func);
void *safe_calloc(size_t nm, size_t sz, int line, char *file, const char *func);
void *safe_realloc(void *ptr, size_t sz, int line, char *file, const char *func);

#undef malloc
#define malloc(x) safe_malloc(x, __LINE__, __FILE__, __func__ )

#undef calloc
#define calloc(n, x) safe_calloc(n, x, __LINE__, __FILE__, __func__ )

#undef realloc
#define realloc(p, x) safe_realloc(p, x, __LINE__, __FILE__, __func__ )

#endif /* INCLUDE_MM_H_ */
