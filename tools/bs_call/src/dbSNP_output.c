/*
 * dbSNP_output.c
 *
 *  Created on: Feb 2, 2020
 *      Author: heath
 */

#define _GNU_SOURCE

#include <sys/types.h>
#include <sys/wait.h>
#include <zlib.h>

#include "dbSNP_idx.h"
#include "dbSNP_utils.h"

static unsigned char dtab2[256] =
{
		33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 0, 0, 0, 0, 0, 133,
		43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 0, 0, 0, 0, 0, 134,
		53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 0, 0, 0, 0, 0, 135,
		63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 0, 0, 0, 0, 0, 136,
		73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 0, 0, 0, 0, 0, 137,
		83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 0, 0, 0, 0, 0, 138,
		93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 0, 0, 0, 0, 0, 139,
		103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 0, 0, 0, 0, 0, 140,
		113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 0, 0, 0, 0, 0, 141,
		123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 0, 0, 0, 0, 0, 142,
};

static int cmp_bin_entries(const void *s1, const void *s2, void *pt) {
	const int ix = *(const int *)s1;
	const int iy = *(const int *)s2;
	const bin *b = pt;
	return (b->entry[ix << 1] & 63) - (b->entry[iy << 1] & 63);
}

#define wrt2bufn(a, sz, buf) { \
	reserve_buffer(buf, sz); \
	memcpy((buf)->mem + (buf)->len, a, sz); \
	(buf)->len += sz; \
}

#define write_buf(f, buf) fwrite((buf)->mem, 1, (buf)->len, f)
#define wrt2buf(a, buf) wrt2bufn(&a, sizeof(a), buf)

static void compress_buf(comp_block_t * const cb, dbsnp_param_t * const par) {
	buffer_t * const buf = &cb->ublock;
	buffer_t * const cbuf = &cb->cblock;
	const size_t size = buf->len;
	if(size > par->max_buf_size) {
		pthread_mutex_lock(&par->param_mut);
		if(size > par->max_buf_size) par->max_buf_size = size;
		pthread_mutex_unlock(&par->param_mut);
	}
	uLong req_size = compressBound((uLong) size);
	resize_buffer(cbuf, req_size);
	uLongf compress_size = cbuf->size;
	int ret = compress((Bytef *)cbuf->mem, &compress_size, (Bytef *)buf->mem, (uLong)size);
	if(ret != 0) {
		fprintf(stderr, "Failed to compress data block\n");
		exit(-3);
	}
	cbuf->len = compress_size;
}

void *compress_thread(void *pt) {
	dbsnp_param_t * const par = pt;
	int idx = 0;
	const int nb = par->n_comp_blocks;
	pthread_mutex_lock(&par->cblock_mut);
	for(;;) {
		for(;;) {
			for(int i = 0; i < nb; i++, idx = (idx + 1) % nb) if(par->cblocks[idx].state == block_uncompressed) break;
			if(par->cblocks[idx].state == block_uncompressed || par->output_finished) break;
			pthread_cond_wait(&par->cblock_cond[1], &par->cblock_mut);
		}
		comp_block_t * const cb = par->cblocks + idx;
		if(cb->state != block_uncompressed) {
			pthread_mutex_unlock(&par->cblock_mut);
			break;
		}
		cb->state = block_active;
		pthread_mutex_unlock(&par->cblock_mut);
		compress_buf(cb, par);
		pthread_mutex_lock(&par->cblock_mut);
		cb->state = block_compressed;
		pthread_cond_signal(&par->cblock_cond[2]);
		idx = (idx + 1) % nb;
	}
	return NULL;
}

void *write_thread(void *pt) {
	dbsnp_param_t * const par = pt;
	int idx = 0;
	const int nb = par->n_comp_blocks;
	FILE *fp = par->outfile;
	pthread_mutex_lock(&par->cblock_mut);
	for(;;) {
		for(;;) {
			if(par->cblocks[idx].state == block_compressed || (par->output_finished && par->cblocks[idx].state == block_empty)) break;
			pthread_cond_wait(&par->cblock_cond[2], &par->cblock_mut);
		}
		comp_block_t * const cb = par->cblocks + idx;
		if(cb->state != block_compressed) {
			pthread_mutex_unlock(&par->cblock_mut);
			break;
		}
		pthread_mutex_unlock(&par->cblock_mut);
		if(cb->flag & 1) *cb->off_ptr = ftell(fp); // First block for contig
		uint64_t sz = cb->cblock.len;
		fwrite(&sz, 1, sizeof(sz), fp);
		write_buf(fp, &cb->cblock);
		if(cb->flag & 2) { // Last block for contig
			// End of contig marker
			uint64_t zero = 0;
			fwrite(&zero, 1, sizeof(zero), fp);
		}
		pthread_mutex_lock(&par->cblock_mut);
		cb->state = block_empty;
		pthread_cond_signal(&par->cblock_cond[0]);
		idx = (idx + 1) % nb;
	}
	return NULL;
}

bool open_outfile(dbsnp_param_t * const par) {
	par->outfile = par->output_file ? fopen(par->output_file, "wb") : stdout;
	if(!par->outfile) {
		fprintf(stderr, "Couldn't open file %s for output: %s\n", par->output_file ? par->output_file : "<stdout>", strerror(errno));
		return true;
	}
	// Skip header block
	fseek(par->outfile, 32, SEEK_SET);
	return false;
}

void finish_output(dbsnp_param_t * const par) {
	FILE * const fp = par->outfile;
	uint64_t off = ftell(fp);
	comp_block_t *cb = par->cblocks;
	buffer_t * const buf = &cb->ublock;
	clear_buffer(buf);
	const uint8_t vers = 2
			;
	const uint8_t zerob = 0;
	wrt2buf(vers, buf);
	wrt2buf(zerob, buf);
	wrt2buf(par->n_prefix, buf);
	uint32_t n_ctgs = HASH_COUNT(par->contigs);
	wrt2buf(n_ctgs, buf);
	if(par->header != NULL) {
		wrt2bufn(par->header, 1 + strlen(par->header), buf);
	} else {
		const char *hd = "track name = dbSNP_index description = \"dbSNP index produced by dbSNP_idx\"";
		wrt2bufn(hd, 1 + strlen(hd), buf);
	}
	for(prefix *p = par->prefixes; p != NULL; p = p->hh.next) {
		wrt2bufn(p->pref, 1 + strlen(p->pref), buf);
	}
	for(contig *ctg = par->contigs; ctg != NULL; ctg = ctg->hh.next) {
		wrt2buf(ctg->min_bin, buf);
		wrt2buf(ctg->max_bin, buf);
		wrt2buf(ctg->offset, buf);
		wrt2bufn(ctg->rname, 1 + strlen(ctg->rname), buf);
	}
	compress_buf(cb, par);
	write_buf(fp, &cb->cblock);
	const uint32_t magic = IDX_MAGIC;
	fwrite(&magic, 1, sizeof(uint32_t), fp);
	fseek(fp, 0, SEEK_SET);
	const uint32_t reserve = 0;
	fwrite(&magic, 1, sizeof(uint32_t), fp);
	fwrite(&reserve, 1, sizeof(uint32_t), fp);
	fwrite(&off, 1, sizeof(uint64_t), fp);
	fwrite(&par->max_buf_size, 1, sizeof(uint64_t), fp);
	uint64_t sz = cb->cblock.len;
	fwrite(&sz, 1, sizeof(sz), fp);
	fclose(fp);

}

static void add_to_compress_queue(buffer_t *buf, dbsnp_param_t * const par, uint8_t flag, uint64_t *off_ptr) {

	const int idx = par->cblock_idx;
	pthread_mutex_lock(&par->cblock_mut);
	while(par->cblocks[idx].state != block_empty) {
		pthread_cond_wait(&par->cblock_cond[0], &par->cblock_mut);
	}
	comp_block_t *cb = par->cblocks + idx;
	swap_buffers(buf, &cb->ublock);
	cb->state = block_uncompressed;
	cb->flag = flag;
	cb->off_ptr = off_ptr;
	pthread_cond_broadcast(&par->cblock_cond[1]);
	pthread_mutex_unlock(&par->cblock_mut);
	clear_buffer(buf);
	par->cblock_idx = (idx + 1) % par->n_comp_blocks;
}

void output_contig(contig *const ctg, dbsnp_param_t *par) {
	int idx[64];
	uint16_t start_ix[64];

	buffer_t *buf = par->output_buf;
	static prefix **pref_list = NULL;
	static int n_prefix = 0;

	pthread_mutex_lock(&par->param_mut);
	if(n_prefix != par->n_prefix) {
		n_prefix = par->n_prefix;
		if(pref_list) free(pref_list);
		pref_list = malloc(sizeof(void *) * n_prefix);
		int ix = 0;
		for(prefix *p = par->prefixes; p != NULL; p = p->hh.next) {
			pref_list[p->ix] = p;
			p->ix = ix++;
		}
	}
	pthread_mutex_unlock(&par->param_mut);
	fprintf(stderr, "Writing output %s\n", ctg->rname);
	int n_items = 0;
	bin *b = ctg->bins;
	uint32_t curr_bin = ctg->min_bin;
	clear_buffer(buf);
	uint8_t fg = 1;
	uint64_t *off_ptr = &ctg->offset;
	uint64_t tn = 0;
	for(uint32_t i = ctg->min_bin; i <= ctg->max_bin; i++, b++) {
		int ne = b->n_entries;
		if(!ne) continue;
		tn += ne;
		// Write out distance from previous bin to new bin
		uint32_t k = i - curr_bin;
		uint8_t x = 0;
		if(k < 64) {
			x = k << 2;
			wrt2buf(x, buf);
		} else if(k < 256) {
			x = 1;
			wrt2buf(x, buf);
			x = k;
			wrt2buf(x, buf);
		} else if(k < 65536) {
			x = 2;
			uint16_t k1 = k;
			wrt2buf(x, buf);
			wrt2buf(k1, buf);
		} else {
			x = 3;
			wrt2buf(x, buf);
			wrt2buf(k, buf);
		}
		curr_bin = i;
		// Sort entries in bin
		for(int j = 0; j < ne; j++) idx[j] = j;
		qsort_r(idx, ne, sizeof(int), cmp_bin_entries, b);
		// Get starting positions for names
		uint16_t x2 = 0;
		for(int j = 0; j < ne; j++) {
			start_ix[j] = x2;
			x2 += b->entry[j << 1] >> 8;
		}
		// Write out entries in bin
		uint8_t terminator = 0;
		for(int j = 0; j < ne; j++) {
			if(j) wrt2buf(terminator, buf);
			int j1 = idx[j];
			if(b->fq_mask & (1l << j1)) {
				terminator = 2;
				par->n_snps_maf_filtered++;
			} else terminator = 0;
			uint16_t z = b->entry[j1 << 1];
			int l = z >> 8;
			uint16_t ix = pref_list[(int)b->entry[(j1 << 1) + 1]]->ix;
			if(ix < 3) z |= ((ix + 1) << 6);
			uint8_t xb = z & 0xff;
			wrt2buf(xb, buf);
			if(ix >= 3) wrt2buf(ix, buf);
			unsigned char *np = b->name_buf + start_ix[j1];
			for(int k = 0; k < l; k++) wrt2bufn(&dtab2[np[k]], 1, buf);
		}
		terminator |= 1;
		wrt2buf(terminator, buf);
		if(++n_items == ITEMS_PER_BLOCK) {
			add_to_compress_queue(buf, par, fg, off_ptr);
			off_ptr = NULL;
			fg = 0;
			n_items = 0;
		}
		free(b->entry);
		free(b->name_buf);
		b->entry = NULL;
		b->name_buf = NULL;
	}
	free(ctg->bins);
	if(n_items > 0) add_to_compress_queue(buf, par, fg | 2, off_ptr);
}

void *output_thread(void *pt) {
	dbsnp_param_t * const par = pt;
	while(true) {
		pthread_mutex_lock(&par->contig_queue_mut);
		while(!par->input_finished && !par->contig_queue) pthread_cond_wait(&par->contig_queue_cond, &par->contig_queue_mut);
		contig *clist = par->contig_queue;
		if(clist) par->contig_queue = clist->next;
		pthread_mutex_unlock(&par->contig_queue_mut);
		if(clist) output_contig(clist, par);
		else if(par->input_finished) break;
	}
	return NULL;
}
