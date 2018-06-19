/*
* compress.c
*
*  Created on: 15 Sep 2016
*      Author: heath
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <locale.h>
#include <pthread.h>
#include <assert.h>
#include <errno.h>

#include "compress.h"
#include "utils.h"

static pthread_mutex_t compress_lock;
static struct compress compress_data;

static void free_compress(void) {
  int i, j;

  if (compress_data.initialized) {
    pthread_mutex_lock(&compress_lock);
    if (compress_data.initialized) {
      for (i = 0; i < COMPRESS_NONE; i++) {
        free(compress_data.compress_suffix[i]);
        for (j = 0; j < 2; j++)
          if (compress_data.comp_path[i][j])
            free(compress_data.comp_path[i][j]);
      }
      compress_data.initialized = false;
    }
    pthread_mutex_unlock(&compress_lock);
  }
}

static void init_compress(void) {
  int i, j;
  char *pnames[][2] = {
      {"bgzip", NULL}, {"gzip", NULL}, {"bzip2", NULL}, {"xz", NULL}, {"compress", NULL}, {NULL, NULL}};
  int compress_type[] = {COMPRESS_GZIP, COMPRESS_GZIP, COMPRESS_BZIP2, COMPRESS_XZ, COMPRESS_COMPRESS, COMPRESS_NONE};
  char *suff[] = {"gz", "bz2", "xz", "Z"};
  char *path;

	if (!compress_data.initialized) {
		pthread_mutex_lock(&compress_lock);
		errno = 0;
		if (!compress_data.initialized) {
			(void)setlocale(LC_ALL, "");
			if (!(path = getenv("PATH")))
			  path = DEFAULT_PATH;
			for (i = 0; i < COMPRESS_NONE; i++) {
				compress_data.compress_suffix[i] = strdup(suff[i]);
				compress_data.comp_path[i][0] = compress_data.comp_path[i][0] = NULL;
			}
			int ix = 0;
			while(pnames[ix][0] != NULL) {
				i = compress_type[ix];
				if(compress_data.comp_path[i][0] == NULL) {
					for (j = 0; j < 2; j++) 
					  compress_data.comp_path[i][j] = pnames[ix][j] ? find_prog(pnames[ix][j], path) : NULL;
				}
				ix++;
			}
			for (i = 0; i < COMPRESS_NONE; i++) if (compress_data.comp_path[i][0] != NULL) break;
			compress_data.default_compress = i;			
			if (atexit(free_compress))
			  fprintf(stderr, "Warning: Unable to register exit function free_compress()\n");
			compress_data.initialized = true;
		}
		errno = 0;
		pthread_mutex_unlock(&compress_lock);
	}
}

struct compress* get_compress_data(void) {
	init_compress();
	return &compress_data;
}
