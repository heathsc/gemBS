/*
 * load_dbSNP.c
 *
 *  Created on: Nov 24, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>
#include <zlib.h>

#include "uthash.h"
#include "dbSNP.h"

static void store_dbsnp_entries(dbsnp_bin_t *bin, int n_entries, int name_buf_sz, uint16_t *entries, uint8_t *name_buf) {
	bin->entries = malloc(sizeof(uint16_t) * n_entries);
	bin->name_buf = malloc((size_t)name_buf_sz);
	bin->n_entries = n_entries;
	uint64_t msk = (uint64_t)0;
	for(int i = 0; i < n_entries; i++) {
		bin->entries[i] = entries[i];
		msk |= ((uint64_t)1 << (entries[i] & 63));
	}
	bin->mask = msk;
	memcpy(bin->name_buf, name_buf, name_buf_sz);
}

dbsnp_header_t *load_dbSNP_header(char * const fname) {
	bool ok = true;
	dbsnp_header_t *hdr = NULL;
	FILE * const file = fopen(fname, "rb");
	if(!file) {
		fprintf(stderr, "Could not open file %s for input: %s\n", fname, strerror(errno));
		return NULL;
	}
	fprintf(stderr,"Loading dbSNP header from %s\n", fname);
	uint32_t td[2];
	size_t sz = fread(td, sizeof(uint32_t), 2, file);
	if(sz != 2 || td[0] != 0xd7278434) {
		fprintf(stderr, "Invalid format\n");
		return NULL;
	}
	void *ucomp_buf = NULL, *comp_buf = NULL;
	uint32_t n_ctgs = 0;
	uint64_t td1[3];
	dbsnp_ctg_t * ctgs = NULL;
	sz = fread(td1, sizeof(uint64_t), 3, file);
	if(sz != 3) ok = false;
	else {
		hdr = calloc(1, sizeof(dbsnp_header_t));
		hdr->fp = file;
		hdr->filename = fname;
		hdr->dbSNP_bufsize = td1[1];
		ucomp_buf = malloc(td1[1]);
		comp_buf = malloc(td1[2]);
		if(fseek(file, td1[0], SEEK_SET)) ok = false;
		else {
			sz = fread(comp_buf, 1, td1[2], file);
			if(sz != td1[2]) ok = false;
			else {
				sz = fread(td, sizeof(uint32_t), 1, file);
				if(sz != 1 || *td != 0xd7278434) ok = false;
			}
		}
		unsigned long size = td1[1];
		if(ok) {
			int ret = uncompress(ucomp_buf, &size, comp_buf, td1[2]);
			if(ret) ok = false;
		}
		if(ok) {
			hdr->n_dbSNP_prefixes = *((uint16_t *)(ucomp_buf + 2));
			n_ctgs = *((uint32_t *)(ucomp_buf + 4));
			ctgs = calloc(n_ctgs, sizeof(dbsnp_ctg_t));
			char *p = ucomp_buf + 8;
			char *p1 = ucomp_buf + size;
			size_t l = strlen(p);
			if(p + 8 >= p1 || strncmp(p, "track ", 6)) ok = false;
			else {
				hdr->dbSNP_header = malloc(l - 5);
				memcpy(hdr->dbSNP_header, p + 6, l - 5);
				hdr->dbSNP_prefix = malloc(sizeof(void *) * hdr->n_dbSNP_prefixes);
				p += l + 1;
			}
			for(int i = 0; ok && i < hdr->n_dbSNP_prefixes && p < p1; i++) {
				l = strlen(p);
				if(p + l >= p1) ok = false;
				else {
					hdr->dbSNP_prefix[i] = malloc(l + 1);
					memcpy(hdr->dbSNP_prefix[i], p, l + 1);
					p += l + 1;
				}
			}
			uint32_t min_bin = 0, max_bin = 0;
			uint64_t offset = 0;
			for(int i = 0; ok && i < n_ctgs && p < p1; i++) {
				if(p + 16 >= p1) ok = false;
				else {
					memcpy(&min_bin, p, sizeof(uint32_t));
					memcpy(&max_bin, p + 4, sizeof(uint32_t));
					memcpy(&offset, p + 8, sizeof(uint64_t));
					if(max_bin < min_bin) {
						ok = false;
					}
					else p += 16;
				}
				if(!ok) break;
				l = strlen(p);
				if(p + l >= p1) ok = false;
				else {
					ctgs[i].min_bin = min_bin;
					ctgs[i].max_bin = max_bin;
					ctgs[i].file_offset = offset;
					ctgs[i].name = malloc(l + 1);
					memcpy(ctgs[i].name, p, l + 1);
					p += l + 1;
				}
			}
		}
	}
	if(ok) {
		for(int i = 0; i < n_ctgs; i++) {
			dbsnp_ctg_t *ctg;
			HASH_FIND(hh, hdr->dbSNP, ctgs[i].name, strlen(ctgs[i].name), ctg);
			if(ctg != NULL) {
				fprintf(stderr,"Error in dbSNP file - duplicate contigs (%s)\n", ctgs[i].name);
				ok = false;
				break;
			}
			HASH_ADD_KEYPTR(hh, hdr->dbSNP, ctgs[i].name, strlen(ctgs[i].name), ctgs + i);
		}
	}
	if(comp_buf) free(comp_buf);
	if(ucomp_buf) free(ucomp_buf);
	if(ok) fprintf(stderr, "dbSNP index loaded OK\n");
	else {
		fprintf(stderr, "dbSNP index loading failed\n");
		free(hdr);
		fclose(file);
		hdr = NULL;
	}
	return hdr;
}

void unload_dbSNP_ctg(dbsnp_ctg_t * const ctg) {
	if(ctg && ctg->bins) {
		dbsnp_bin_t *bin = ctg->bins;
		for(uint32_t bn = ctg->min_bin; bn <= ctg->max_bin; bn++, bin++) {
			if(bin->n_entries) {
				free(bin->entries);
				free(bin->name_buf);
			}
		}
		free(ctg->bins);
		ctg->bins = NULL;
	}
}

bool load_dbSNP_ctg(const dbsnp_header_t * const hdr, dbsnp_ctg_t * const ctg) {
	bool ok = true;

	static uint8_t db_tab[] = {
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x10, 0x11, 0x12, 0x13, 0x14,
			0x15, 0x16, 0x17, 0x18, 0x19, 0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x30,
			0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46,
			0x47, 0x48, 0x49, 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x60, 0x61, 0x62,
			0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78,
			0x79, 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x90, 0x91, 0x92, 0x93, 0x94,
			0x95, 0x96, 0x97, 0x98, 0x99, 0x0f, 0x1f, 0x2f, 0x3f, 0x4f, 0x5f, 0x6f, 0x7f, 0x8f, 0x9f, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff
	};

	fprintf(stderr,"Loading dbSNP entries for %s\n", ctg->name);
	FILE * const file = hdr->fp;

	void * const ucomp_buf = malloc(hdr->dbSNP_bufsize);
	size_t comp_buf_size = 1 + hdr->dbSNP_bufsize * .75;
	void *comp_buf = malloc(comp_buf_size);
	if(fseek(file, ctg->file_offset, SEEK_SET)) ok = false;
	uint16_t *entries = malloc(sizeof(uint16_t) * 64);
	uint8_t *name_buf = malloc(sizeof(uint8_t) * 256 * 64);
	int n_snps = 0, n_bins = 0;
	ctg->bins = calloc((ctg->max_bin - ctg->min_bin + 1), sizeof(dbsnp_bin_t));
	dbsnp_bin_t *bins = ctg->bins;
	uint32_t curr_bin = ctg->min_bin;
	while(ok) {
		uint64_t sz;
		size_t k = fread(&sz, sizeof(uint64_t), 1, file);
		unsigned long size;
		if(k != 1) ok = false;
		else if(sz == 0) break;
		else {
			if(comp_buf_size < sz) {
				comp_buf_size = sz * 1.1;
				comp_buf = realloc(comp_buf, comp_buf_size);
			}
			k = fread(comp_buf, 1, sz, file);
			if(k != sz) ok = false;
			else {
				size = hdr->dbSNP_bufsize;
				int ret = uncompress(ucomp_buf, &size, comp_buf, sz);
				if(ret) ok = false;
			}
		}
		if(!ok) break;
		uint8_t *bp = ucomp_buf;
		uint8_t *bp_end = ucomp_buf + size;
		int n_entries = 0, name_buf_ptr = 0;
		bool end_of_bin = false;
		int prev_ix = -1;
		while(ok && bp < bp_end) {
			if(!n_entries) {
				uint32_t bin_inc = 0;
				uint8_t x = *bp++;
				switch(x & 3) {
				case 0:
					bin_inc = x >> 2;
					break;
				case 1:
					if(bp < bp_end) {
						x = *bp++;
						bin_inc = x;
					} else ok = false;
					break;
				case 2:
					if(bp + 1 < bp_end) {
						uint16_t k;
						memcpy(&k, bp, 2);
						bp += 2;
						bin_inc = k;

					} else ok = false;
					break;
				case 3:
					if(bp + 3 < bp_end) {
						uint32_t k;
						memcpy(&k, bp, 4);
						bp += 4;
						bin_inc = k;
					} else ok = false;
					break;
				}
				if(!ok) break;
				curr_bin += bin_inc;
				if(curr_bin > ctg->max_bin || bp >= bp_end) break;
				bins += bin_inc;
			}
			uint8_t x = *bp++;
			int prefix_ix = x >> 6;
			int sl;
			if(!prefix_ix) {
				if(bp + 2 < bp_end) {
					name_buf[name_buf_ptr++] = *bp++;
					name_buf[name_buf_ptr++] = *bp++;
					sl = 2;
				} else ok = false;
			} else sl = 0;
			if(!ok) break;
			if((x & 63) <= prev_ix || prefix_ix > hdr->n_dbSNP_prefixes) ok = false;
			else {
				prev_ix = x & 63;
				int k = name_buf_ptr;
				while((*bp) > 1 && sl++ < 256 && bp < bp_end) {
					name_buf[name_buf_ptr++] = db_tab[(int)(*bp++)];
				}
				k = name_buf_ptr - k;
				if(*bp > 1) ok = false;
				else {
					if(*bp++ == 1) end_of_bin = true;
					if(n_entries == 64) ok = false;
					else entries[n_entries++] = (k << 8) | (uint16_t)x;
				}
			}
			if(!ok) break;
			if(end_of_bin) {
				store_dbsnp_entries(bins, n_entries, name_buf_ptr, entries, name_buf);
				n_bins++;
				n_snps += n_entries;
				n_entries = 0;
				name_buf_ptr = 0;
				prev_ix = -1;
				end_of_bin = false;
			}
		}
	}
	fprintf(stderr, "ctg loaded: %s, n_snps = %d\n", ok ? "OK" : "BAD", n_snps);
	free(ucomp_buf);
	free(comp_buf);
	free(entries);
	free(name_buf);
	return ok;
}

bool dbSNP_lookup_name(const dbsnp_header_t *const hdr, const dbsnp_ctg_t * ctg, char * const rs, size_t * const rs_len, const uint32_t x) {
	static char dtab[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 0, 0, 0, 0, 0, 0 };
	bool found = false;
	rs[0] = 0;
	if(ctg != NULL) {
		int bn = x >> 6;
		if(bn >= ctg->min_bin && bn <= ctg->max_bin) {
			dbsnp_bin_t *b = ctg->bins + bn - ctg->min_bin;
			int ix = x & 63;
			uint64_t mk = (uint64_t)1 << ix;
			if(b->mask & mk) {
				uint64_t mk1 = b->mask & (mk - (uint64_t)1);
				int i = 0, j = 0;
				while(mk1) {
					if(mk1 & (uint64_t)1) {
						uint16_t en = b->entries[i++];
						j += en >> 8;
						if(!((en >> 6) & 3)) j += 2;
					}
					mk1 >>= 1;
				}
				char *tp = rs;
				int prefix_id = (b->entries[i] >> 6) & 3;
				unsigned char *tp1 = b->name_buf + j;
				if((prefix_id--) == 0) {
					prefix_id = (tp1[0] << 8) | tp1[1];
					tp1+=2;
				}
				char * db_prefix = hdr->dbSNP_prefix[prefix_id];
				j = 0;
				while(db_prefix[j]) *tp++ = db_prefix[j++];
				j = b->entries[i] >> 8;
				for(int k = 0; k < j; k++) {
					unsigned char z = *tp1++;
					*tp++ = dtab[z >> 4];
					*tp++ = dtab[z & 15];
				}
				*tp = 0;
				if(rs_len) *rs_len = tp - rs;
				found = true;
			}
		}
	}
	return found;
}
