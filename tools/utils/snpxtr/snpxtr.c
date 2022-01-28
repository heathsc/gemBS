#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <pthread.h>
#include <time.h>
#include <sys/wait.h>
#include "utils.h"
#include "snpxtr.h"

#include "htslib/hfile.h"
#include "htslib/bgzf.h"

void error(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
	exit(-1);
}

int main(int argc, char **argv) {

	sargs_t args;
	init_params(&args);
	handle_command_line(argc, argv, &args);
	int ns = bcf_hdr_nsamples(args.hdr);
	assert(ns > 0);
	if(args.snplistname) read_snplist(&args);
	bool error = false;
	if(args.dbSNPfilename) {
		args.dbSNP_hdr = load_dbSNP_header(args.dbSNPfilename);
		if(!args.dbSNP_hdr) error = true;
	}
	if(!error) {
		process_input(&args);
		fprintf(stderr,"snpxtr: finished\n");
	}
	bcf_sr_destroy(args.sr);
	return 0;
}


