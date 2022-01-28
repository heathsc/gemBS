/*
 * dbSNP_command_line.c
 *
 *  Created on: Feb 1, 2020
 *      Author: heath
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <stdbool.h>

#include "dbSNP_idx.h"

static const char *usage(void) {
	return
			"\n"
			"About: Create index of dbSNP positions for bs_call.\n"
			"Usage: dbSNP_idx [-o outfile] [infile] [infile] ...\n"
			"Options:\n"
			"   -o, --outfile <FILE>                Output file for index (default, stdout)\n"
			"   -d, --desc                         Description of dataset\n"
			"   -t, --type <AUTO, BED, JSON, VCF>  Input file type (default, BED)\n"
			"   -u, --unsorted-file <FILE>          Input file that has unsorted records from multiple contigs\n"
			"   -c, --chrom-alias <FILE>           Chromosome name alias file\n"
			"   -S, --sorted                       Assume input files are sorted by contigs (unless specified with -u command)\n"
			"   -@, --threads <n>                  Extra threads\n"
			"   -m, --maf-limit <x>                Flag SNPs with maf >= x so that bs_call will output these sites even if\n"
			"                                      they are AA or TT reference homozygous.  This option only works with JSON input files\n"
			"   -s, --selected <FILE>              File with list of SNPs (one per line) to be flagged as for --maf-limit\n"
			"\n";
}

static struct option loptions[] = {
		{"outfile",required_argument,0,'o'},
		{"desc",required_argument,0,'d'},
		{"type",required_argument,0,'t'},
		{"unsorted-file",required_argument,0,'u'},
		{"chrom-alias",required_argument,0,'c'},
		{"threads",required_argument,0,'@'},
		{"sorted",no_argument,0,'S'},
		{"maf",no_argument,0,'m'},
		{"selected",no_argument,0,'s'},
		{0,0,0,0}
};

void add_file(dbsnp_param_t * const par, char * const name, const bool sorted) {
	file_t *us, **base;
	if(sorted) {
		HASH_FIND_STR(par->unsorted, name, us);
		if(us) return;
		base = &par->sorted;
	} else base = &par->unsorted;
	HASH_FIND_STR((*base), name, us);
	if(!us) {
		us= malloc(sizeof(file_t));
		us->processed = false;
		us->name = name;
		us->sorted = sorted;
		us->read = false;
		HASH_ADD_KEYPTR(hh, (*base), us->name, strlen(us->name), us);
	}
}

void handle_command_line(int argc, char *argv[], dbsnp_param_t * const par) {
	int c;
	bool sorted = false;
	while ((c = getopt_long(argc, argv, "o:u:c:t:m:s:d:@:Sh?",loptions,NULL)) >= 0) {
		switch (c) {
		case 'o':
			par->output_file = optarg;
			break;
		case 'd':
			par->header = optarg;
			break;
		case 'c':
			par->chrom_alias_file = optarg;
			break;
		case 's':
			par->select_file = optarg;
			break;
		case 't':
			if(!strcasecmp(optarg, "BED")) par->input_type = dbsnp_bed;
			else if(!strcasecmp(optarg, "JSON")) par->input_type = dbsnp_json;
			else if(!strcasecmp(optarg, "VCF")) par->input_type = dbsnp_vcf;
			else if(!strcasecmp(optarg, "AUTO")) par->input_type = dbsnp_auto;
			break;
		case 'u':
			add_file(par, optarg, false);
			break;
		case 'S':
			sorted = true;
			break;
		case 'm': {
			double z = atof(optarg);
			if(z < 0.0 || z > 0.5) fprintf(stderr, "--maf-limit must be between 0 and 0.5\n");
			else par->maf_limit = z;
		}
		break;
		case '@':
			par->threads = atoi(optarg);
			if(par->threads < 0) par->threads = 0;
			break;
		case 'h':
		case '?':
		default:
			fputs(usage(), stdout);
			exit(0);
			break;
		}
	}
	if(!par->threads) par->threads = sysconf(_SC_NPROCESSORS_ONLN) - 1;
	for(int i = optind; i <= argc; i++) {
		if(i == argc) {
			if(i == optind) add_file(par, "-", sorted);
		} else add_file(par, argv[i], sorted);
	}
	if(par->select_file) read_select_file(par);
}
