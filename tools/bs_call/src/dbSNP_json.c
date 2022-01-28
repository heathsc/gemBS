/*
 * dbSNP_json.c
 *
 *  Created on: Feb 10, 2020
 *      Author: heath
 */

#include "dbSNP_idx.h"
#include <ctype.h>

#include "jsmn.h"

typedef enum {
	refsnp_id, primary_snapshot_data, placements_with_allele, is_ptlp, alleles,
	allele, spdi, inserted_sequence, deleted_sequence, position, seq_id,
	allele_annotations, frequency, study_name, allele_count, total_count, observation
} keys_t;

#define IN_PRIMARY_SNAPSHOT_DATA 1
#define IN_PLACEMENTS_WITH_ALLELE 2
#define IN_ALLELES 4
#define IN_ALLELE 8
#define IN_SPDI 0x10
#define IN_ALLELE_ANNOTATIONS 0x20
#define IN_FREQUENCY 0x40
#define IN_OBSERVATION 0x80
#define IS_PTLP 0x100
#define HAS_POSITION 0x200
#define FREQ_ALLELES_OK 0x400
#define VALID_SNP 0x800
#define SEEN_FREQ_INSERTED_SEQUENCE 0x1000
#define SEEN_FREQ_DELETED_SEQUENCE 0x2000
#define SEEN_ALLELE_COUNT 0x4000
#define SEEN_TOTAL_COUNT 0x8000
#define STUDY_NAME_OK 0x10000

#define FREQ_ALLELE_FLAGS (SEEN_FREQ_INSERTED_SEQUENCE | SEEN_FREQ_DELETED_SEQUENCE)
#define FREQ_FLAGS (FREQ_ALLELES_OK | SEEN_ALLELE_COUNT | SEEN_TOTAL_COUNT | STUDY_NAME_OK)

typedef struct {
	keys_t val;
	char *key;
	UT_hash_handle hh;
} key_hash_t;

typedef struct {
	snp_t *snp;
	key_hash_t *h;
	uint32_t allele_count;
	uint32_t total_count;
	uint32_t a;
	uint32_t b;
	uint32_t mask;
	char inserted_sequence;
	char deleted_sequence;
	char alleles[2];
} snp_extr_t;

typedef struct {
	key_hash_t *keys;
	jsmn_parser *p;
	size_t tcount;
	jsmntok_t *tok;
} jsmn_work_t;

static void init_jsmn_work(jsmn_work_t * const work) {
	char *keywords[] = {
			"refsnp_id",
			"primary_snapshot_data",
			"placements_with_allele",
			"is_ptlp",
			"alleles",
			"allele",
			"spdi",
			"inserted_sequence",
			"deleted_sequence",
			"position",
			"seq_id",
			"allele_annotations",
			"frequency",
			"study_name",
			"allele_count",
			"total_count",
			"observation",
			NULL
			};

	work->keys = NULL;
	int ct = 0;
	keys_t vals[] = {
			refsnp_id, primary_snapshot_data, placements_with_allele, is_ptlp, alleles,
			allele, spdi, inserted_sequence, deleted_sequence, position, seq_id,
			allele_annotations, frequency, study_name, allele_count, total_count, observation
	};
	while(keywords[ct]) {
		key_hash_t *key;
		const size_t sz = strlen(keywords[ct]);
		HASH_FIND(hh, work->keys, keywords[ct], sz, key);
		if(!key) {
			key = malloc(sizeof(key_hash_t));
			key->key = keywords[ct];
			key->val = vals[ct];
			HASH_ADD_KEYPTR(hh, work->keys, key->key, sz, key);
		}
		ct++;
	}
	work->tcount = 1024;
	work->tok = malloc (sizeof(jsmntok_t) * work->tcount);
	work->p = malloc(sizeof(jsmn_parser));
}

static int handle_json(char * const js, jsmntok_t * const tok, int count, snp_extr_t * const snp, const int level) {
	int i, j;
	key_hash_t *kh;

	if(!count || (snp->mask & VALID_SNP)) return 0;
	switch(tok->type) {
	case JSMN_PRIMITIVE:
	case JSMN_STRING:
		return 1;
		break;
	case JSMN_OBJECT:
		j = 0;
		for(i = 0; i < tok->size && !(snp->mask & VALID_SNP); i++) {
			jsmntok_t *ntok = tok + j + 1;
			assert(ntok->type == JSMN_STRING);
			HASH_FIND(hh, snp->h, js + ntok->start, ntok->end - ntok->start, kh);
			bool child_processed = false;
			if(kh && ntok->size > 0) {
				switch(kh->val) {
				case refsnp_id:
					if(!level && ntok[1].type == JSMN_STRING) {
						snp->snp->name = js + ntok[1].start;
						snp->snp->name_len = ntok[1].end - ntok[1].start;
						j += 2;
						child_processed = true;
					}
					break;
				case primary_snapshot_data:
					if(!level && ntok[1].type == JSMN_OBJECT) {
						snp->mask |= IN_PRIMARY_SNAPSHOT_DATA;
						j++;
						j += handle_json(js, tok + j + 1, count - j, snp, level + 1);
						snp->mask &= ~IN_PRIMARY_SNAPSHOT_DATA;
						if(snp->inserted_sequence && snp->deleted_sequence && snp->inserted_sequence != snp->deleted_sequence && (snp->mask & HAS_POSITION)) {
							snp->mask |= VALID_SNP;
							if(snp->total_count > 0) {
								double z = (double)snp->allele_count / (double)snp->total_count;
								if(z > 0.5) z = 1.0 - z;
								snp->snp->maf = z;
							}
						}
						child_processed = true;
					}
					break;
				case placements_with_allele:
					if(level == 1 && ntok[1].type == JSMN_ARRAY && (snp->mask & IN_PRIMARY_SNAPSHOT_DATA)) {
						snp->mask |= IN_PLACEMENTS_WITH_ALLELE;
						j++;
						j += handle_json(js, tok + j + 1, count - j, snp, level + 1);
						snp->mask &= ~IN_PLACEMENTS_WITH_ALLELE; // Clear all allele flags
						child_processed = true;
					}
					break;
				case allele_annotations:
					if(level == 1 && ntok[1].type == JSMN_ARRAY && (snp->mask & IN_PRIMARY_SNAPSHOT_DATA)) {
						snp->mask |= IN_ALLELE_ANNOTATIONS;
						j++;
						j += handle_json(js, tok + j + 1, count - j, snp, level + 1);
						snp->mask &= ~IN_ALLELE_ANNOTATIONS; // Clear all allele flags
						child_processed = true;
					}
					break;
				case is_ptlp:
					if(level == 3 && ntok[1].type == JSMN_PRIMITIVE && (snp->mask & IN_PLACEMENTS_WITH_ALLELE)) {
						if(js[ntok[1].start] == 't') snp->mask |= IS_PTLP;
						else snp->mask &= ~IS_PTLP;
						j += 2;
						child_processed = true;
					}
					break;
				case frequency:
					if(level == 3 && ntok[1].type == JSMN_ARRAY && (snp->mask & IN_ALLELE_ANNOTATIONS)) {
						snp->mask = (snp->mask & ~FREQ_FLAGS) | IN_FREQUENCY;
						j += 2;
						for(int i1 = 0; i1 < ntok[1].size; i1++) {
							j += handle_json(js, tok + j + 1, count - j, snp, level + 2);
							if((snp->mask & FREQ_FLAGS) == FREQ_FLAGS && snp->a <= snp->b) {
								snp->allele_count += snp->a;
								snp->total_count += snp->b;
							}
							snp->mask &= ~FREQ_FLAGS;
						}
						snp->mask &= ~IN_FREQUENCY;
						child_processed = true;
					}
					break;
				case alleles:
					if(level == 3 && ntok[1].type == JSMN_ARRAY && ((snp->mask & (IN_PLACEMENTS_WITH_ALLELE | IS_PTLP)) == (IN_PLACEMENTS_WITH_ALLELE | IS_PTLP))) {
						snp->mask |= IN_ALLELES;
						j++;
						j += handle_json(js, tok + j + 1, count - j, snp, level + 1);
						snp->mask &= ~IN_ALLELES;
						child_processed = true;
					}
					break;
				case allele:
					if(level == 5 && ntok[1].type == JSMN_OBJECT && (snp->mask & IN_ALLELES)) {
						snp->mask |= IN_ALLELE;
						snp->inserted_sequence = snp->deleted_sequence = 0;
						uint32_t old_pos = snp->snp->pos;
						uint16_t snp_mask = snp->mask & HAS_POSITION;
						j++;
						j += handle_json(js, tok + j + 1, count - j, snp, level + 1);
						if(snp->inserted_sequence && snp->deleted_sequence && snp->inserted_sequence != snp->deleted_sequence && (snp->mask & HAS_POSITION)) {
							//							snp->mask |= VALID_SNP;
							snp->alleles[0] = snp->deleted_sequence;
							snp->alleles[1] = snp->inserted_sequence;
						} else {
							snp->snp->pos = old_pos;
							snp->mask = (snp->mask & ~(HAS_POSITION)) | snp_mask;
						}
						snp->mask &= ~IN_ALLELE; // Clear all allele flags
						child_processed = true;
					}
					break;
				case observation:
					if(level == 5 && ntok[1].type == JSMN_OBJECT && (snp->mask & IN_FREQUENCY)) {
						snp->mask = (snp->mask & ~FREQ_ALLELE_FLAGS) | (IN_OBSERVATION | FREQ_ALLELES_OK);
						j++;
						j += handle_json(js, tok + j + 1, count - j, snp, level + 1);
						if((snp->mask & (FREQ_ALLELE_FLAGS | FREQ_ALLELES_OK)) != (FREQ_ALLELE_FLAGS | FREQ_ALLELES_OK)) snp->mask &= ~(FREQ_ALLELE_FLAGS | FREQ_ALLELES_OK);
						else snp->mask &= ~FREQ_ALLELE_FLAGS;
						child_processed = true;
					}
					break;
				case spdi:
					if(level == 6 && ntok[1].type == JSMN_OBJECT && (snp->mask & IN_ALLELE)) {
						snp->mask |= IN_SPDI;
						j++;
						j += handle_json(js, tok + j + 1, count - j, snp, level + 1);
						snp->mask &= ~IN_SPDI; // Clear all allele flags
						child_processed = true;
					}
					break;
				case position:
					if(level == 7 && ntok[1].type == JSMN_PRIMITIVE && (snp->mask & IN_SPDI)) {
						char *tp = js + ntok[1].start;
						if(isdigit((int)*tp)) {
							char *tp1;
							snp->snp->pos = strtoul(tp, &tp1, 10);
							if(tp1 == js + ntok[1].end) {
								snp->mask |= HAS_POSITION;
							}
						}
						j += 2;
						child_processed = true;
					}
					break;
				case seq_id:
					if(level == 7 && ntok[1].type == JSMN_STRING && (snp->mask & IN_SPDI)) {
						snp->snp->cname = js + ntok[1].start;
						snp->snp->cname_len = ntok[1].end - ntok[1].start;
						j += 2;
						child_processed = true;
					}
					break;
				case inserted_sequence:
					if(level == 7 && ntok[1].type == JSMN_STRING && (snp->mask & IN_SPDI)) {
						if(ntok[1].end - ntok[1].start == 1) {
							snp->inserted_sequence = js[ntok[1].start];
						}
						j += 2;
						child_processed = true;
					} else if(level == 6 && ntok[1].type == JSMN_STRING && (snp->mask & IN_OBSERVATION)) {
						snp->mask |= SEEN_FREQ_INSERTED_SEQUENCE;
						if(ntok[1].end - ntok[1].start != 1 || js[ntok[1].start] != snp->inserted_sequence) snp->mask &= ~FREQ_ALLELES_OK;
						j += 2;
						child_processed = true;
					}
					break;
				case deleted_sequence:
					if(level == 7 && ntok[1].type == JSMN_STRING && (snp->mask & IN_SPDI)) {
						if(ntok[1].end - ntok[1].start == 1) {
							snp->deleted_sequence = js[ntok[1].start];
						}
						j += 2;
						child_processed = true;
					} else if(level == 6 && ntok[1].type == JSMN_STRING && (snp->mask & IN_OBSERVATION)) {
						snp->mask |= SEEN_FREQ_DELETED_SEQUENCE;
						if(ntok[1].end - ntok[1].start != 1 || js[ntok[1].start] != snp->deleted_sequence) snp->mask &= ~FREQ_ALLELES_OK;
						j += 2;
						child_processed = true;
					}
					break;
				case study_name:
					if(level == 5 && ntok[1].type == JSMN_STRING && (snp->mask & IN_FREQUENCY)) {
						snp->mask |= STUDY_NAME_OK;
						j += 2;
						child_processed = true;
					}
					break;
				case allele_count:
					if(level == 5 && ntok[1].type == JSMN_PRIMITIVE && (snp->mask & IN_FREQUENCY)) {
						char *tp = js + ntok[1].start;
						if(isdigit((int)*tp)) {
							char *tp1;
							snp->a = strtoul(tp, &tp1, 10);
							if(tp1 == js + ntok[1].end) {
								snp->mask |= SEEN_ALLELE_COUNT;
							}
						}
						j += 2;
						child_processed = true;
					}
					break;
				case total_count:
					if(level == 5 && ntok[1].type == JSMN_PRIMITIVE && (snp->mask & IN_FREQUENCY)) {
						char *tp = js + ntok[1].start;
						if(isdigit((int)*tp)) {
							char *tp1;
							snp->b = strtoul(tp, &tp1, 10);
							if(tp1 == js + ntok[1].end) {
								snp->mask |= SEEN_TOTAL_COUNT;
							}
						}
						j += 2;
						child_processed = true;
					}
					break;
				default:
					break;
				}
			}
			if(!child_processed) {
				j++;
				if(ntok->size > 0) {
					j += handle_json(js, tok + j + 1, count - j, snp, level + 1);
				}
			}
		}
		return j + 1;
		break;
	case JSMN_ARRAY:
		j = 0;
		for(i = 0; i < tok->size && !(snp->mask & VALID_SNP); i++) {
			j += handle_json(js, tok + j + 1, count - j, snp, level + 1);
		}
		return j + 1;
	default:
		return 0;
		break;
	}
}

void parse_json_line(char * const buf, const ssize_t l, jsmn_work_t * const work, snp_t * const snp, dbsnp_param_t * const par) {
	if(!work->keys) init_jsmn_work(work);
	int r;
	jsmn_init(work->p);
	snp_extr_t s;
	do {
		r = jsmn_parse(work->p, buf, l, work->tok, work->tcount);
		if(r < 0) {
			if(r == JSMN_ERROR_NOMEM) work->tcount = work->tcount << 1;
			work->tok = realloc(work->tok, sizeof(jsmntok_t) * work->tcount);
		} else break;
	} while(r < 0);
	if(r > 0) {
		memset(&s, 0, sizeof(s));
		s.snp = snp;
		s.h = work->keys;
		handle_json(buf, work->tok, r, &s, 0);
		if((s.mask & VALID_SNP) && snp->name && snp->cname) {
//			fprintf(stdout, "rs%.*s\t%.*s\t%u\t%c\t%c\t%g\n", snp->name_len, snp->name, snp->cname_len, snp->cname, snp->pos, s.alleles[0], s.alleles[1], snp->maf);
			snp->ok = true;
		}
	}
}
