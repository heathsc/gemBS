/*
 * bbi_structs.h
 *
 *  Created on: Jan 10, 2020
 *      Author: heath
 */

#ifndef BBI_STRUCTS_H_
#define BBI_STRUCTS_H_

#include "bbi_defs.h"

// Main header of bbi file
typedef struct {
	uint32_t magic;
	uint16_t version;
	uint16_t zoomLevels;
	uint64_t chromosomeTreeOffset;
	uint64_t fullDataOffset;
	uint64_t fullIndexOffset;
	uint16_t fieldCount;
	uint16_t definedFieldCount;
	uint64_t autoSqlOffset;
	uint64_t totalSummaryOffset;
	uint32_t uncompressBufSize;
	uint64_t extensionOffset;
} bbi_header_t;

typedef struct {
	uint32_t start;
	uint32_t end;
	uint64_t offset;
} bbi_block_t;

typedef enum { cblock_empty = 0, cblock_uncompressed, cblock_active, cblock_compressed } bbi_cblock_state_t;

typedef struct {
	int ctg_id;
	int block_idx;
	kstring_t *buf_p;
	int ix;
	bbi_cblock_state_t state;
} bbi_cblock_t;

typedef struct {
	uint32_t end_base;
	uint32_t count;
} bb_zrec_t;

typedef struct {
	uint32_t end_base;
	uint32_t count;
	float x;
	float xsq;
	float min;
	float max;
} bw_zrec_t;

typedef struct {
	union {
		bb_zrec_t *bb_rec;
		bw_zrec_t *bw_rec;
	};
	int size;
	int ix;
} bbi_zblock_t;

typedef struct {
	bw_zrec_t *bw_rec;
	int size;
	int ix;
} bw_zblock_t;

typedef struct {
	uint32_t start;
	float val;
} bw_rec_t;

typedef struct {
	bbi_block_t bbuf;
	bbi_block_t *blocks;
	bbi_zblock_t zblock[ZOOM_LEVELS - 1];
	bw_rec_t bw_rec[BW_ITEMS_PER_SLOT];
	int n_items;
	int block_idx;
	int block_sz;
} bbi_data_t;

// Stored data to allow generation of zoom levels for bigBed and bigWig files
// We store data on two bases in each byte of base_type;
// bits 4-7: base 1, bits 0-3: base 2
//
// bits 0-1: bedmethyl_type
// bit 2: strand (0 == top, 1 == bottom)
// bit 3: non-zero methylation (1 == yes)
// bits 4-5: bedmethyl_type
// bit 6: strand (0 == top, 1 == bottom)
// bit 7: non-zero methylation (1 == yes)
//
// There is one zoom_data_t structure per contig
//

typedef struct {
	uint8_t *base_type;
	float *val;
	int val_size;
	int val_ix;
	uint64_t len; // Chromosome length
} zoom_dt_t;

typedef struct {
	zoom_dt_t zoom_data;
	bbi_data_t bbi_data[5];
} bbi_ctg_data_t;

typedef struct {
	kstring_t *buffer;
	void *comp_buf;
	uint32_t zoom_scales[ZOOM_LEVELS];
	size_t comp_buf_size;
	uint32_t n_rec;
	uint64_t index_offset;
	uint64_t max_buf_size;
	uint64_t zoom_data_offset[ZOOM_LEVELS];
	uint64_t zoom_index_offset[ZOOM_LEVELS];
	uint32_t first_base;
	uint32_t last_base;
	uint32_t first_ctg;
	uint32_t last_ctg;
	uint32_t n_cts;
	uint32_t res_size[ZOOM_LEVELS];
	uint32_t res_end[ZOOM_LEVELS];
	uint64_t total_bases;
	double min_x;
	double max_x;
	double sum_x;
	double sum_xsq;
	bool first_time;
} bbi_global_data_t;


#endif /* BBI_STRUCTS_H_ */
