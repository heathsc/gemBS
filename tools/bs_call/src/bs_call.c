/*
 * bs_call.c
 *
 *  Created on: 30 Sep 2014
 *      Author: heath
 */

#define BS_CALL "bs_call"

#include <stdio.h>
#include <sys/wait.h>
#include "gem_tools.h"

#include "bs_call.h"

int main(int argc, char *argv[]) {
	gt_status err = 0;
	// GT error handler
	gt_handle_error_signals();
	sr_param param;

	init_param(&param);

	// Parsing command-line options
	err = parse_arguments(argc, argv, &param);
	if (err == GT_STATUS_OK) {
		err = bs_call_process(&param);
	}
	fprintf(stderr, "Finished processing\n");
	if(param.work.json_file != NULL) {
		fprintf(stderr, "Writing out statistics\n");
		output_stats(&param);
		fclose(param.work.json_file);
	}
	while(waitpid(-1, NULL, 0) > 0);
	fprintf(stderr, "bs_call finishing\n");
	return err == GT_STATUS_OK ? 0 : -1;
}
