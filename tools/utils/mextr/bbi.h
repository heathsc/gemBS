/*
 * bbi.h
 *
 *  Created on: Jan 7, 2020
 *      Author: heath
 */



#ifndef BBI_H_
#define BBI_H_

#include "bbi_structs.h"

#define bbi_write(stream,src) fwrite(&src, sizeof(src), 1, stream)

typedef struct {
	uint64_t n; // Number of leaf level nodes
	uint64_t n_items; // Number of items
	uint32_t block_size;
	int depth; // Depth of tree (including leaves)
	int wsize; // Size of width array
	int *width; // Number of nodes at each level
	void **start; // Starting position for each node (size is width + 1 with last element giving total number of nodes)
} tree_t;

// Information about the contig B+ tree
typedef struct {
	char **names;
	uint64_t *len;
	uint32_t key_len;
	tree_t tree;
} contig_tree_t;

typedef struct {
	uint32_t start_ctg;
	uint32_t end_ctg;
	uint32_t start_base;
	uint32_t end_base;
	int start_idx; // Starting index for node (or data block) at lower level
} r_node_t;

// Information about the R trees (main index + zoom indices)
typedef struct {
	uint32_t ctg;
	bbi_block_t * block;
} r_tree_block_t;

typedef struct {
	int nctgs;
	uint64_t end_offset;
	r_tree_block_t *blocks;
	tree_t tree;
} r_tree_t;

typedef struct {
	args_t * args;
	int ix;
	uint32_t nrec;
} bbi_thr_info_t;

void init_bbi_header(args_t * const args, const bool bigbed);
void write_bbi_header(FILE * const fp, bbi_header_t * const header, bbi_global_data_t * const bdata);
void finish_bb_block(args_t * const args, const int ctg_id, const int ix);
void finish_bw_block(args_t * const args, const int ctg_id, const int ix);
void finish_bbi_blocks(args_t * const args, const int ctg_id);
void finish_bb_data_file(args_t * const args, const int ix);
void finish_bw_data_file(args_t * const args, const int ix);
void *bbi_write_thread(void *p);
void *bbi_compress_thread(void *p);
void *handle_bedmethyl_thread(void *p);
void *bbi_create_index(void *p);
void init_cblocks(args_t * const args, const int nb);
void clear_cblocks(args_t * const args);
void destroy_cblocks(args_t * const args);

#endif /* BBI_H_ */
