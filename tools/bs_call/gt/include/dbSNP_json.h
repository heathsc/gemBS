/*
 * dbSNP_json.h
 *
 *  Created on: Feb 10, 2020
 *      Author: heath
 */

#ifndef INCLUDE_DBSNP_JSON_H_
#define INCLUDE_DBSNP_JSON_H_

#include "dbSNP_idx.h"

typedef struct {
	void *keys;
	void *p;
	void *tok;
	size_t tcount;
} jsmn_work_t;

void parse_json_line(char * const buf, const ssize_t l, jsmn_work_t *work, snp_t * const snp, dbsnp_param_t * const par);

#endif /* INCLUDE_DBSNP_JSON_H_ */
