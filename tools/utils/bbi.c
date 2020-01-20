/*
 * bbi.c
 *
 *  Created on: Jan 7, 2020
 *      Author: heath
 *
 *  Routines to write bigBed and bigWig files from the bedmethyl and wig files
 *  generated for the standard ENCODE analysis pipeline.
 *
 *  This is not and is not intended to be a general purpose library for writing
 *  bigBed and bigWig files, and is instead a specialized set of routines that take
 *  advantage of the peculiarities of the data (all ranges of size 1, most cytosines
 *  not methylated etc.) to increase speed and memory efficiency.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <pthread.h>
#include <zlib.h>

#include "htslib/khash_str2int.h"
#include "mextr.h"
#include "bbi.h"

static const char *autosql_desc =
		"table BisulfiteSeq\n"
		"\"BED9+5 scores for bisulfite-seq data\"\n"
		"\t(\n"
		"\tstring\tchrom;\t\"Reference chromosome or scaffold\"\n"
		"\tuint\tchromStart;\t\"Start position in chromosome\"\n"
		"\tuint\tchromEnd;\t\"End position in chromosome\"\n"
		"\tstring\tname;\t\"Name of item\"\n"
		"\tuint\tscore;\t\"Score from 0-1000.  Capped number of reads\"\n"
		"\tchar[1]\tstrand;\t\"+ or - or . for unknown\"\n"
		"\tuint\tthickStart;\t\"Start of where display should be thick (start codon)\"\n"
		"\tuint\tthickEnd;\t\"End of where display should be thick (stop codon)\"\n"
		"\tuint\treserved;\t\"Color value R,G,B\"\n"
		"\tuint\treadCount;\t\"Number of reads or coverage\"\n"
		"\tuint\tpercentMeth;\t\"Percentage of reads that show methylation at this position in the genome\"\n"
		"\tstring\trefContext;\t\"Reference context on strand (2 bases for CpG, 3 bases for CHG, CHH)\"\n"
		"\tstring\tcalledContext;\t\"Called context on strand (2 bases for CpG, 3 bases for CHG, CHH)\"\n"
		"\tuint\tgenotypeQual;\t\"Phred score for genotype call\"\n"
		"\t)\n";

void write_bbi_header(FILE * const fp, bbi_header_t * const header, bbi_global_data_t * const bdata) {
	assert(fp);
	fseek(fp, 0, SEEK_SET);
	bbi_write(fp, header->magic);
	bbi_write(fp, header->version);
	bbi_write(fp, header->zoomLevels);
	bbi_write(fp, header->chromosomeTreeOffset);
	bbi_write(fp, header->fullDataOffset);
	bbi_write(fp, header->fullIndexOffset);
	bbi_write(fp, header->fieldCount);
	bbi_write(fp, header->definedFieldCount);
	bbi_write(fp, header->autoSqlOffset);
	bbi_write(fp, header->totalSummaryOffset);
	bbi_write(fp, header->uncompressBufSize);
	bbi_write(fp, header->extensionOffset);
	// zoomHeaders
	uint32_t res = 0;
	for(int i = 0; i < ZOOM_LEVELS; i++) {
		bbi_write(fp, bdata->zoom_scales[i]);
		bbi_write(fp, res);
		bbi_write(fp, bdata->zoom_data_offset[i]);
		bbi_write(fp, bdata->zoom_index_offset[i]);
	}
	// Add autoSql string
	if(header->autoSqlOffset) {
		fseek(fp, header->autoSqlOffset, SEEK_SET);
		fwrite(autosql_desc, 1, 1 + strlen(autosql_desc), fp);
	}
	// And total summary
	bbi_write(fp, bdata->total_bases);
	bbi_write(fp, bdata->min_x);
	bbi_write(fp, bdata->max_x);
	bbi_write(fp, bdata->sum_x);
	bbi_write(fp, bdata->sum_xsq);
	// And extended header
	fseek(fp, header->extensionOffset, SEEK_SET);
	uint16_t ext_size = 64;
	bbi_write(fp, ext_size);
	for(int i = 0; i < 62; i++) fputc(0, fp);
}

static int ctg_id_lookup(const tree_t * tree, const int i, const int level) {
	const int j = ((const int * const)tree->start[level])[i];
	return level > 0 ? ctg_id_lookup(tree, j, level - 1) : j;
}

void write_contig_tree_level(FILE *fp, const contig_tree_t *ctree, const int level) {
	uint64_t pos = ftell(fp);
	const tree_t * tree = &ctree->tree;
	char * const keybuf = malloc(ctree->key_len);
	const int n_nodes = tree->width[level]; // Number of nodes at this level
	const uint64_t item_size = ctree->key_len + 8; // item size is the same for leaf and non-leaf nodes
	const int * const start_ix = tree->start[level];
	if(level > 0) {
		// Non-leaf nodes
		const int n_nodes1 = tree->width[level - 1]; // Number of nodes at next (lower) level
		uint64_t off = pos + 4 * n_nodes + item_size * n_nodes1; // Offset of first node at next (lower) level
		const int * const start_ix1 = tree->start[level - 1];
		uint16_t zero = 0;
		for(int i = 0; i < n_nodes; i++) {
			uint16_t n_items = start_ix[i+1] - start_ix[i];
			bbi_write(fp, zero);
			bbi_write(fp, n_items);
			for(int j = start_ix[i]; j < start_ix[i + 1]; j++) {
				// Look up chromosome name corresponding to item position
				uint32_t id = ctg_id_lookup(tree, j, level - 1);
				strncpy(keybuf, ctree->names[(int)id], ctree->key_len);
				fwrite(keybuf, 1, ctree->key_len, fp);
				bbi_write(fp, off);
				off += 4 + item_size * (start_ix1[j + 1] - start_ix1[j]);
			}
		}
	} else {
		// Leaf nodes
		for(int i = 0; i < n_nodes; i++) {
			uint8_t isLeaf = 1;
			uint8_t res = 0;
			uint16_t n_items = start_ix[i + 1] - start_ix[i];
			bbi_write(fp, isLeaf);
			bbi_write(fp, res);
			bbi_write(fp, n_items);
			for(int j = start_ix[i]; j < start_ix[i + 1]; j++) {
				strncpy(keybuf, ctree->names[j], ctree->key_len);
				uint32_t id = j;
				uint32_t csize = ctree->len[id];
				fwrite(keybuf, 1, ctree->key_len, fp);
				bbi_write(fp, id);
				bbi_write(fp, csize);
			}
		}
	}
	free(keybuf);
}

void write_contig_tree(FILE * const fp, bbi_header_t * const header, const contig_tree_t * const ctree) {
	assert(fp);
	fseek(fp, header->chromosomeTreeOffset, SEEK_SET);
	// Write chromosome B+ tree header
	uint32_t magic = 0x78CA8C91;
	uint32_t valSize = 8;
	uint64_t reserved = 0;
	bbi_write(fp, magic);
	bbi_write(fp, ctree->tree.block_size);
	bbi_write(fp, ctree->key_len);
	bbi_write(fp, valSize);
	bbi_write(fp, ctree->tree.n_items);
	bbi_write(fp, reserved);
	// Write contig tree levels
	for(int i = ctree->tree.depth - 1; i >= 0; i--) write_contig_tree_level(fp, ctree, i);
}

void write_r_tree_level(FILE *fp, const r_tree_t *rtree, const int level) {
	uint64_t pos = ftell(fp);
	const tree_t * tree = &rtree->tree;
	const int n_nodes = tree->width[level]; // Number of nodes at this level
	const uint64_t leaf_item_size = 32;
	const uint64_t non_leaf_item_size = 24;
	const r_node_t * const rn = tree->start[level];
	if(level > 0) {
		// Non-leaf nodes
		const int n_nodes1 = tree->width[level - 1]; // Number of nodes at next (lower) level
		uint64_t off = pos + 4 * n_nodes + non_leaf_item_size * n_nodes1; // Offset of first node at next (lower) level
		const r_node_t * const rn1 = tree->start[level - 1];
		uint16_t zero = 0;
		int item_size1 = level > 1 ? non_leaf_item_size : leaf_item_size;
		for(int i = 0; i < n_nodes; i++) {
			uint16_t n_items = rn[i + 1].start_idx - rn[i].start_idx;
			bbi_write(fp, zero);
			bbi_write(fp, n_items);
			for(int j = rn[i].start_idx; j < rn[i + 1].start_idx; j++) {
				bbi_write(fp, rn1[j].start_ctg);
				bbi_write(fp, rn1[j].start_base);
				bbi_write(fp, rn1[j].end_ctg);
				bbi_write(fp, rn1[j].end_base);
				bbi_write(fp, off);
				off += 4 + item_size1 * (rn1[j + 1].start_idx - rn1[j].start_idx);
			}
		}
	} else {
		// Leaf nodes
		r_tree_block_t * const rbp = rtree->blocks;
		for(int i = 0; i < n_nodes; i++) {
			uint8_t isLeaf = 1;
			uint8_t res = 0;
			uint16_t n_items = rn[i + 1].start_idx - rn[i].start_idx;
			bbi_write(fp, isLeaf);
			bbi_write(fp, res);
			bbi_write(fp, n_items);
			for(int j = rn[i].start_idx; j < rn[i + 1].start_idx; j++) {
				bbi_write(fp, rbp[j].ctg);
				bbi_write(fp, rbp[j].block->start);
				bbi_write(fp, rbp[j].ctg);
				bbi_write(fp, rbp[j].block->end);
				bbi_write(fp, rbp[j].block->offset);
				uint64_t size = (j + 1 < tree->n ? rbp[j + 1].block->offset : rtree->end_offset) - rbp[j].block->offset;
				bbi_write(fp, size);
			}
		}
	}
}

void write_r_tree(args_t const * args, const int ix, r_tree_t * const rtree) {
	FILE * const fp = ix < 3 ? args->bigbedfiles[ix] : args -> bigwigfiles[ix - 3];
	const uint32_t magic = 0x2468ACE0;
	bbi_write(fp, magic);
	bbi_write(fp, rtree->tree.block_size);
	bbi_write(fp, rtree->tree.n_items);
	r_node_t * const root_node = rtree->tree.start[rtree->tree.depth - 1];
	bbi_write(fp, root_node->start_ctg);
	bbi_write(fp, root_node->start_base);
	bbi_write(fp, root_node->end_ctg);
	bbi_write(fp, root_node->end_base);
	bbi_write(fp, args->bb_global[ix].index_offset);
	const uint32_t items_per_slot = ix < 3 ? ITEMS_PER_SLOT : BW_ITEMS_PER_SLOT;
	const uint32_t res = 0;
	bbi_write(fp, items_per_slot);
	bbi_write(fp, res);
	// Write R tree levels
	for(int i = rtree->tree.depth - 1; i >= 0; i--) write_r_tree_level(fp, rtree, i);
}

void calc_start_ix(tree_t * const tree, const int level, const bool rtree, const r_tree_block_t * const rbp) {
	// For higher level nodes we balance node sizes across the tree, but for the level 1 nodes
	// (one above the leaves) they all have to have the same number of entries (block_size) apart
	// from the last one which normally has less.  This is to allow the reader to quickly go from
	// chromosome ID to the key
	int w = tree->width[level];
	int *start_ix = NULL;
	r_node_t *rn = NULL, *rn1 = NULL;
	if(rtree) {
		rn = tree->start[level];
		if(level > 0) rn1 = tree->start[level - 1];
	}
	else start_ix = tree->start[level];
	int w1 = level > 0 ? tree->width[level - 1] : tree->n;
	if(rtree || level > 0) {
		int k = w1 / w;
		int o1 = w * k - w1;
		int o2 = w * (k + 1) - w1;
		int d = 0;
		int off = 0;
		int sz;
		for(int i = 0; i < w; i++) {
			if(abs(d + o1) < abs(d + o2)) {
				d += o1;
				sz = k;
			} else {
				d += o2;
				sz = k + 1;
			}
			if(rtree) {
				rn[i].start_idx = off;
				if(rn1) {
					rn[i].start_ctg = rn1[off].start_ctg;
					rn[i].start_base = rn1[off].start_base;
					rn[i].end_ctg = rn1[off + sz - 1].end_ctg;
					rn[i].end_base = rn1[off + sz -1].end_base;
				} else {
					rn[i].start_ctg = rbp[off].ctg;
					rn[i].start_base = rbp[off].block->start;
					rn[i].end_ctg = rbp[off + sz - 1].ctg;
					rn[i].end_base = rbp[off + sz - 1].block->end;
				}
			} else start_ix[i] = off;
			off += sz;
		}
		assert(off == w1);
	} else for(int i = 0; i < w; i++) start_ix[i] = i * tree->block_size;
	if(rtree) rn[w].start_idx = w1;
	else start_ix[w] = w1;
}

void set_tree_widths(tree_t * const tree, uint32_t const block_size, size_t item_size, const bool rtree, const r_tree_block_t * const rbp) {
	uint64_t n1 = tree->n;
	while(true) {
		hts_resize(int, tree->depth + 1, &tree->wsize, &tree->width, 0);
		n1 = (n1 + BLOCK_SIZE - 1) / BLOCK_SIZE;
		tree->width[tree->depth++] = n1;
		if(n1 <= 1) break;
	}
	tree->block_size = tree->depth > 1 ? block_size : tree->n;
	tree->start = calloc(tree->depth, sizeof(void *));
	for(int i = 0; i < tree->depth; i++) {
		tree->start[i] = malloc(item_size * (tree->width[i] + 1));
		calc_start_ix(tree, i, rtree, rbp);
	}
}

contig_tree_t *init_contig_tree(args_t * const args) {
	contig_tree_t * const ctree = calloc(1, sizeof(contig_tree_t));
	// Pick up contig names from the selected regions
	// These are already in alphabetical order and filtered so that
	// only contigs contained in the BCF/VCF header are present
	bcf_sr_regions_t * const reg = args->sr->regions;
	const int n = ctree->tree.n = ctree->tree.n_items = reg->nseqs;
	assert(n > 0);
	ctree->len = malloc(n * sizeof(uint64_t));
	ctree->names = reg->seq_names;
	int max = 0;
	for(int i = 0; i < n; i++) {
		const int j = bcf_hdr_name2id(args->hdr, ctree->names[i]);
		assert(j >= 0);
		const bcf_idpair_t * const idp = args->hdr->id[BCF_DT_CTG] + j;
		ctree->len[i] = idp->val->info[0];
		int s = (int)strlen(ctree->names[i]);
		if(s > max) max = s;
	}
	ctree->key_len = max;
	tree_t * const tree = &ctree->tree;
	set_tree_widths(tree, BLOCK_SIZE, sizeof(int), false, NULL);
	return ctree;
}

r_tree_t *init_r_tree(args_t * const args, const int ix, const uint32_t n_items, uint64_t data_end){
	r_tree_t * const rtree = calloc(1, sizeof(r_tree_t));
	bcf_sr_regions_t * const reg = args->sr->regions;
	rtree->nctgs = reg->nseqs;
	rtree->end_offset = data_end;
	tree_t * const tree = &rtree->tree;
	uint64_t nb = 0;
	for(int i = 0; i < rtree->nctgs; i++) nb += args->ctg_data[i].bbi_data[ix].block_idx;
	tree->n = nb;
	tree->n_items = n_items;
	rtree->blocks = malloc(sizeof(r_tree_block_t) * nb);
	nb = 0;
	r_tree_block_t *rbp = rtree->blocks;
	for(int i = 0; i < rtree->nctgs; i++) {
		bbi_data_t * const bbd = args->ctg_data[i].bbi_data + ix;
		for(int j = 0; j < bbd->block_idx; j++, rbp++) {
			rbp->block = bbd->blocks + j;
			rbp->ctg = i;
		}
	}
	set_tree_widths(tree, BLOCK_SIZE, sizeof(r_node_t), true, rtree->blocks);
	return rtree;
}

void destroy_tree(tree_t * const tree) {
	if(tree->width) free(tree->width);
	if(tree->start) {
		for(int i = 0; i < tree->depth; i++) free(tree->start[i]);
		free(tree->start);
	}
}

void destroy_contig_tree(contig_tree_t * const ctree) {
	destroy_tree(&ctree->tree);
	free(ctree);
}

void destroy_r_tree_and_data(args_t const * args, const int ix, r_tree_t * const rtree) {
	destroy_tree(&rtree->tree);
	if(rtree->blocks) free(rtree->blocks);
	for(int i = 0; i < rtree->nctgs; i++) {
		bbi_data_t * const bbd = args->ctg_data[i].bbi_data + ix;
		if(bbd->blocks) free(bbd->blocks);
		bbd->block_idx = bbd->block_sz = bbd->n_items = 0;
		bbd->blocks = NULL;
		memset(&bbd->bbuf, 0, sizeof(bbd->bbuf));
	}
	free(rtree);

}

void _init_bbi_header(bbi_header_t * const header, const bool bigbed) {
	memset(header, 0, sizeof(bbi_header_t));
	header->magic = bigbed ? 0x8789F2EB : 0x888FFC26;
	header->version = 4;
	header->zoomLevels = ZOOM_LEVELS;
	header->totalSummaryOffset = BBI_HEADER_SIZE + header->zoomLevels * ZOOM_HEADER_SIZE;
	if(bigbed) {
		header->fieldCount = 14;
		header->definedFieldCount = 0;
		header->autoSqlOffset = header->totalSummaryOffset;
		header->totalSummaryOffset += strlen(autosql_desc) + 1;
	}
	header->extensionOffset = header->totalSummaryOffset + TOTAL_SUMMARY_SIZE;
	header->chromosomeTreeOffset = header->extensionOffset + EXT_HEADER_SIZE;
}

void init_bbi_header(args_t * const args, const bool bigbed) {
	const int ix = bigbed ? 0 : 1;
	args->bbi_hdr[ix] = malloc(sizeof(bbi_header_t));
	_init_bbi_header(args->bbi_hdr[ix], bigbed);
	contig_tree_t *ctree = init_contig_tree(args);
	args->ctg_data = calloc(ctree->tree.n, sizeof(bbi_ctg_data_t));
	bbi_ctg_data_t * bdata = args->ctg_data;
	for(int i = 0; i < ctree->tree.n; i++, bdata++) {
		bdata->zoom_data.base_type = calloc(1, (ctree->len[i] + 1) >> 1);
		bdata->zoom_data.len = ctree->len[i];
	}
	if(bigbed) {
		for(int i = 0; i < 3; i++) write_contig_tree(args->bigbedfiles[i], args->bbi_hdr[0], ctree);
		args->bbi_hdr[0]->fullDataOffset = ftell(args->bigbedfiles[0]);
		for(int i = 0; i < 3; i++) fseek(args->bigbedfiles[i], args->bbi_hdr[0]->fullDataOffset + 4, SEEK_SET);
	} else {
		const int j = args->strand_specific ? 2 : 1;
		for(int i = 0; i < j; i++) write_contig_tree(args->bigwigfiles[i], args->bbi_hdr[1], ctree);
		args->bbi_hdr[1]->fullDataOffset = ftell(args->bigwigfiles[0]);
		for(int i = 0; i < j; i++) fseek(args->bigwigfiles[i], args->bbi_hdr[1]->fullDataOffset + 4, SEEK_SET);
	}
	destroy_contig_tree(ctree);
}

void *bbi_write_thread(void *p) {
	uint32_t idx = 0;
	args_t * const args = p;
	cblock_buffer_t * const cbuf = &args->cblock_buf;
	const int nb = cbuf->n_cblocks;
	for(;;) {
		pthread_mutex_lock(&cbuf->mut);
		for(;;) {
			if((cbuf->cblocks[idx].state == cblock_compressed) || (cbuf->end_of_input && cbuf->cblocks[idx].state == cblock_empty)) break;
			pthread_cond_wait(&cbuf->cond[2], &cbuf->mut);
		}
		bbi_cblock_t * const cb = cbuf->cblocks + idx;
		if(cb->state != cblock_compressed) {
			pthread_mutex_unlock(&cbuf->mut);
			break;
		}
		const int ix = cb->ix;
		FILE * const fp = ix < 3 ? args->bigbedfiles[ix] : args->bigwigfiles[ix - 3];
		args->ctg_data[cb->ctg_id].bbi_data[ix].blocks[cb->block_idx].offset = ftell(fp);
		pthread_mutex_unlock(&cbuf->mut);
		kstring_t * buf = cb->buf_p;
		fwrite(ks_str(buf), 1, ks_len(buf), fp);
		cb->state = cblock_empty;
		pthread_cond_signal(&cbuf->cond[0]);
		idx = (idx + 1) % nb;
	}
	pthread_cond_signal(&cbuf->cond[0]);
	return NULL;
}

void *bbi_compress_thread(void *p) {
	uint32_t idx = 0;
	args_t * const args = p;
	void *comp_buf = NULL;
	size_t comp_buf_size = 0;
	cblock_buffer_t * const cbuf = &args->cblock_buf;
	const int nb = cbuf->n_cblocks;
	for(;;) {
		pthread_mutex_lock(&cbuf->mut);
		for(;;) {
			for(int i = 0; i < nb; i++, idx = (idx + 1) % nb) if(cbuf->cblocks[idx].state == cblock_uncompressed) break;
			if(cbuf->cblocks[idx].state == cblock_uncompressed || cbuf->end_of_input) break;
			pthread_cond_wait(&cbuf->cond[1], &cbuf->mut);
		}
		bbi_cblock_t * const cb = cbuf->cblocks + idx;
		if(cb->state != cblock_uncompressed) {
			pthread_mutex_unlock(&cbuf->mut);
			break;
		}
		cb->state = cblock_active;
		pthread_mutex_unlock(&cbuf->mut);
		kstring_t * buf = cb->buf_p;
		const int ix = cb->ix;
		uLong req_size = compressBound((uLong) ks_len(buf));
		if(req_size > comp_buf_size) {
			comp_buf_size = req_size * 1.2;
			comp_buf = realloc(comp_buf, comp_buf_size);
		}
		if(ks_len(buf) > args->bb_global[ix].max_buf_size) args->bb_global[ix].max_buf_size = ks_len(buf);
		uLongf compress_size = comp_buf_size;
		int ret = compress((Bytef *)comp_buf, &compress_size, (Bytef *)ks_str(buf), (uLong)ks_len(buf));
		if(ret != 0) error("Failed to compress data block\n");
		ks_resize(buf, (size_t)compress_size);
		memcpy(buf->s, comp_buf, (size_t)compress_size);
		buf->l = (size_t)compress_size;
		cb->state = cblock_compressed;
		pthread_cond_signal(&cbuf->cond[2]);
		idx = (idx + 1) % nb;
	}
	if(comp_buf) free(comp_buf);
	pthread_cond_signal(&cbuf->cond[2]);
	return NULL;
}

void init_cblocks(args_t * const args, const int nb) {
	args->cblock_buf.n_cblocks = nb;
	args->cblock_buf.cblocks = calloc(nb, sizeof(bbi_cblock_t));
	for(int i = 0; i < args->cblock_buf.n_cblocks; i++) {
		args->cblock_buf.cblocks[i].buf_p = malloc(sizeof(kstring_t));
		ks_initialize(args->cblock_buf.cblocks[i].buf_p);
	}
}

void clear_cblocks(args_t * const args) {
	for(int i = 0; i < args->cblock_buf.n_cblocks; i++) {
		args->cblock_buf.cblocks[i].state = cblock_empty;
		ks_clear(args->cblock_buf.cblocks[i].buf_p);
	}
	args->cblock_buf.end_of_input = false;
	args->cblock_buf.pos = 0;
}

void destroy_cblocks(args_t * const args) {
	if(args->cblock_buf.n_cblocks) {
		for(int i = 0; i < args->cblock_buf.n_cblocks; i++) ks_free(args->cblock_buf.cblocks[i].buf_p);
		free(args->cblock_buf.cblocks);
		args->cblock_buf.cblocks = NULL;
		args->cblock_buf.n_cblocks = 0;
	}
}

void finish_bb_block(args_t * const args, const int ctg_id, const int ix) {
	bbi_global_data_t * const gdata = args->bb_global + ix;
	kstring_t *buf = gdata->buffer;
	if(ks_len(buf) > 0) {
		// Get cblock to store block for compression
		cblock_buffer_t * const cbuf = &args->cblock_buf;
		const int pos = cbuf->pos;
		pthread_mutex_lock(&cbuf->mut);
		while(cbuf->cblocks[pos].state != cblock_empty) {
			pthread_cond_wait(&cbuf->cond[0], &cbuf->mut);
		}
		bbi_data_t * const bdata = args->ctg_data[ctg_id].bbi_data + ix;
		hts_resize(bbi_block_t, bdata->block_idx + 1, &bdata->block_sz, &bdata->blocks, 0);
		bbi_block_t * const bl = bdata->blocks + (bdata->block_idx++);
		memcpy(bl, &bdata->bbuf, sizeof(bdata->bbuf));
		bbi_cblock_t * const bp = cbuf->cblocks + pos;
		bp->ix = ix;
		gdata->buffer = bp->buf_p;
		bp->buf_p = buf;
		bp->ctg_id = ctg_id;
		bp->block_idx = bdata->block_idx - 1;
		bp->state = cblock_uncompressed;
		pthread_mutex_unlock(&cbuf->mut);
		pthread_cond_signal(&cbuf->cond[1]);
		cbuf->pos = (pos + 1) % cbuf->n_cblocks;
		if(gdata->first_time) {
			gdata->first_ctg = ctg_id;
			gdata->first_base = bdata->bbuf.start;
			gdata->first_time = false;
		}
		gdata->last_ctg = ctg_id;
		gdata->last_base = bdata->bbuf.end;
		ks_clear(gdata->buffer);
		memset(&bdata->bbuf, 0, sizeof(bdata->bbuf));
		bdata->n_items = 0;
	}
}

void finish_bw_block(args_t * const args, const int ctg_id, const int ix) {
	bbi_data_t * const bdata = args->ctg_data[ctg_id].bbi_data + 3 + ix;
	if(bdata->n_items > 0) {
		// Write bigWig section to buffer
		bbi_global_data_t * const gdata = args->bb_global + ix + 3;
		gdata->n_rec++;
		bw_rec_t * bw_rec = bdata->bw_rec;
		kstring_t *buf = gdata->buffer;
		uint32_t dat[5] = {ctg_id, bw_rec[0].start, bw_rec[bdata->n_items - 1].start + 1, 0, 1};
		uint8_t dat1[2] = {2, 0};
		uint16_t itemCount = bdata->n_items;
		kputsn_((char *)dat, sizeof(uint32_t) * 5, ks_clear(buf));
		kputsn_((char *)dat1, 2, buf);
		kputsn_((char *)&itemCount, 2, buf);
		for(int i = 0; i < bdata->n_items; i++, bw_rec++) {
			kputsn_((char *)&bw_rec->start, 4, buf);
			kputsn_((char *)&bw_rec->val, 4, buf);
		}
		// Get cblock to store block for compression
		cblock_buffer_t * const cbuf = &args->cblock_buf;
		const int pos = cbuf->pos;
		pthread_mutex_lock(&cbuf->mut);
		while(cbuf->cblocks[pos].state != cblock_empty) {
			pthread_cond_wait(&cbuf->cond[0], &cbuf->mut);
		}
		bbi_data_t * const bdata = args->ctg_data[ctg_id].bbi_data + 3 + ix;
		hts_resize(bbi_block_t, bdata->block_idx + 1, &bdata->block_sz, &bdata->blocks, 0);
		bbi_block_t * const bl = bdata->blocks + (bdata->block_idx++);
		bl->start = dat[1];
		bl->end = dat[2];
		bbi_cblock_t * const bp = cbuf->cblocks + pos;
		bp->ix = ix + 3;
		gdata->buffer = bp->buf_p;
		bp->buf_p = buf;
		bp->ctg_id = ctg_id;
		bp->block_idx = bdata->block_idx - 1;
		bp->state = cblock_uncompressed;
		pthread_mutex_unlock(&cbuf->mut);
		pthread_cond_signal(&cbuf->cond[1]);
		cbuf->pos = (pos + 1) % cbuf->n_cblocks;
		if(gdata->first_time) {
			gdata->first_ctg = ctg_id;
			gdata->first_base = dat[1];
			gdata->first_time = false;
		}
		gdata->last_ctg = ctg_id;
		gdata->last_base = dat[2];
		ks_clear(gdata->buffer);
		memset(&bdata->bbuf, 0, sizeof(bdata->bbuf));
		bdata->n_items = 0;
	}
}

void finish_bbi_blocks(args_t * const args, const int ctg_id) {
	for(int i = 0; i < 3; i++) finish_bb_block(args, ctg_id, i);
	for(int i = 0; i < 2; i++) finish_bw_block(args, ctg_id, i);
}

void finish_bb_data_file(args_t * const args, const int ix) {
	FILE *fp = args->bigbedfiles[ix];
	args->bb_global[ix].index_offset = ftell(fp);
}

void finish_bw_data_file(args_t * const args, const int ix) {
	FILE *fp = args->bigwigfiles[ix];
	args->bb_global[ix + 3].index_offset = ftell(fp);
}

// Create and write main index
void *bbi_create_index(void *p) {
	bbi_thr_info_t * const bi = p;
	args_t * const args = bi->args;
	const int ix = bi->ix;
	// Create and write main index
	uint64_t pos = ftell(ix < 3 ? args->bigbedfiles[ix] : args->bigwigfiles[ix - 3]);
	r_tree_t * const rtree = init_r_tree(args, ix, bi->nrec, pos);
	write_r_tree(args, ix, rtree);
	destroy_r_tree_and_data(args, ix, rtree);
	return NULL;
}
