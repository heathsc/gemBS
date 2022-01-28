/*
 * init_params.c
 *
 *  Created on: Dec 26, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <stdbool.h>
#include <pthread.h>
#include "utils.h"
#include "snpxtr.h"

#include "htslib/hfile.h"

void init_params(sargs_t *const args) {
	memset(args, 0, sizeof(sargs_t));
	ks_initialize(&args->out_string);
}
