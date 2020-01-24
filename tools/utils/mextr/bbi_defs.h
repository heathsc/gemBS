/*
 * bbi_defs.h
 *
 *  Created on: Jan 15, 2020
 *      Author: heath
 */

#ifndef BBI_DEFS_H_
#define BBI_DEFS_H_

#define BLOCK_SIZE 256
#define ITEMS_PER_SLOT 512
#define BW_ITEMS_PER_SLOT 1024

#define INITIAL_REDUCTION 10
#define BW_INITIAL_REDUCTION 40
#define ZOOM_RES_INCREMENT 4

#define BBI_HEADER_SIZE 64
#define EXT_HEADER_SIZE 64
#define ZOOM_HEADER_SIZE 24
#define TOTAL_SUMMARY_SIZE 40

#define ZOOM_LEVELS 10

#endif /* BBI_DEFS_H_ */
