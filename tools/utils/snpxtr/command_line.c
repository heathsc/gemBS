/*
 * command_line.c
 *
 *  Created on: Dec 26, 2019
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

#include "htslib/hfile.h"
#include "htslib/khash_str2int.h"

#include "utils.h"
#include "snpxtr.h"

// These are copied from htslib:synced_bcf_reader.c as the definitions are not visible
// in the standard library and we need them to allow sorting of regions

typedef struct {
	hts_pos_t start, end;
} region1_t;

struct bcf_sr_region_t {
	region1_t *regs;
	int nregs, mregs, creg;
} region_t;

const char *usage(void) {
	return
			"\n"
			"About: Extract SNPs from VCF/BCF file.\n"
			"Usage: snpxtr [file] [regions]\n"
			"Options:\n"
			"   -o, --output            Output file (default = stdout)\n"
			"   -s, --snps              File with list of SNPs to be selected (default, select all sites with PASS)\n"
			"   -D, --dbsnp             dbSNP index file (used to add external ids if not present in input file\n"
			"   -r, --regions           restrict to comma separated list of regions\n"
			"   -R, --regions-file       restrict to regions listed in file\n"
			"   -@, --threads           Extra threads\n"
			"   -z, --bgzip             Compress output with bgzip\n"
			"   -m, --md5               Calculate md5 digest for output file (if not stdout)\n"
			"   -x, --tabix             Generate tabix (tbx) index for compressed output file\n"
			"\n";
}

static struct option loptions[] = {
		{"output",required_argument,0,'o'},
		{"snps",required_argument,0,'s'},
		{"dbsnp",required_argument,0,'D'},
		{"regions",required_argument,0,'r'},
		{"regions-file",required_argument,0,'R'},
		{"regions",required_argument,0,'r'},
		{"threads",required_argument,0,'@'},
		{"bgzip",no_argument,0,'z'},
		{"md5",no_argument,0,'m'},
		{"tabix",no_argument,0,'x'},
		{"help",no_argument,0,'h'},
		{0,0,0,0}
};


void handle_command_line(int argc, char *argv[], sargs_t * const args) {
	int c;
	bool regions_file = false;
	char *regions_list = NULL;
	while ((c = getopt_long(argc, argv, "o:s:D:r:R:@:zmxh?",loptions,NULL)) >= 0) {
		switch (c) {
		case 'o':
			args->outfilename = optarg;
			break;
		case 's':
			args->snplistname = optarg;
			break;
		case 'D':
			args->dbSNPfilename = optarg;
			break;
		case 'x':
			args->tabix = true;
			break;
		case 'm':
			args->md5 = true;
			break;
		case 'R':
			regions_file = true;
			// fall through
		case 'r':
			regions_list = optarg;
			break;
		case '@':
			args->threads = atoi(optarg);
			if(args->threads < 0) args->threads = 0;
			break;
		case 'z':
			args->compress = true;
			break;
		case 'h':
		case '?':
		default: error(usage()); break;
		}
	}
	char *fname = NULL;
	if(optind == argc) error(usage());
	else fname = argv[optind];
	args->sr = bcf_sr_init();
	bcf_sr_set_threads(args->sr, args->threads);
	// Process region arguments if present
	if(regions_list) {
		if(bcf_sr_set_regions(args->sr, regions_list, regions_file) < 0) error("Failed to parse the regions: %s\n", regions_list);
	} else if(optind + 1 < argc) {
		kstring_t tmp = {0, 0, 0};
		kputs(argv[optind + 1], &tmp);
		for(int k = optind + 2; k < argc; k++) {
			kputc(',', &tmp);
			kputs(argv[k], &tmp);
		}
		if(bcf_sr_set_regions(args->sr, tmp.s, 0) < 0) error("Failed to parse the regions: %s\n", tmp.s);
		free(tmp.s);
	}
	if(args->threads > 0) bcf_sr_set_threads(args->sr, args->threads);
	if(!bcf_sr_add_reader(args->sr, fname))
		error("failed to read from %s: %s\n", fname, bcf_sr_strerror(args->sr->errnum));
	args->hdr = args->sr->readers[0].header;
	args->pass_idx = args->gt_idx = -1;
	for(int i = 0; i < args->hdr->n[BCF_DT_ID]; i++) {
		if(args->pass_idx < 0 && !strcasecmp("PASS", args->hdr->id[BCF_DT_ID][i].key)) args->pass_idx = i;
		if(args->gt_idx < 0 && !strcasecmp("GT", args->hdr->id[BCF_DT_ID][i].key)) args->gt_idx = i;
	}
	int ns = bcf_hdr_nsamples(args->hdr);
	assert(ns > 0);
	args->gt = malloc(sizeof(int) * ns);
	bcf_sr_regions_t * const reg = args->sr->regions;
	int nctgs_vcf = args->hdr->n[BCF_DT_CTG];
	args->cumul_len = malloc(nctgs_vcf * sizeof(uint64_t));
	if(reg) {
		for(int i = 0; i < nctgs_vcf; i++) {
			const bcf_idpair_t * const idp = args->hdr->id[BCF_DT_CTG] + i;
			int k;
			uint64_t len = 0;
			int ret = khash_str2int_get(args->sr->regions->seq_hash, idp->key, &k);
			if(ret >= 0) len = reg->regs[k].regs->end + 1;
			args->cumul_len[i] = i ? args->cumul_len[i - 1] + len : len;
		}
		int nctgs = reg->nseqs;
	} else {
		for(int i = 0; i < nctgs_vcf; i++) {
			const bcf_idpair_t * const idp = args->hdr->id[BCF_DT_CTG] + i;
			uint64_t len = idp->val->info[0];
			args->cumul_len[i] = i ? args->cumul_len[i - 1] + len : len;
		}
	}
	if(!args->outfilename) args->outfilename = "-";
	args->outfile = open_ofile(&args->outfilename, args->compress, args);
}
