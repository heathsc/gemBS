/*
 * read_utils.c
 *
 *  Created on: Nov 26, 2019
 *      Author: heath
 */

#include <stdio.h>

#include "gem_tools.h"
#include "bs_call.h"

void trim_read(gt_vector * const sqread, int left_trim, int right_trim) {
	if(sqread != NULL) {
		uint32_t rl = gt_vector_get_used(sqread);
		if(rl > 0) {
			uint8_t *const sp = gt_vector_get_mem(sqread, uint8_t);
			for (int k1 = 0; k1 < left_trim && k1 < rl; k1++) {
				sp[k1] = (sp[k1] & 3) | ((FLT_QUAL) << 2);
			}
			for (int k1 = 0; k1 < right_trim && k1 < rl; k1++)	{
				sp[rl - k1 - 1] = (sp[k1] & 3) | ((FLT_QUAL) << 2);
			}
		}
	}
}
