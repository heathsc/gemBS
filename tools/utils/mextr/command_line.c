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
#include "mextr.h"

// These are copied from htslib:synced_bcf_reader.c as the definitions are not visible
// in the standard library and we need them to allow sorting of regions

typedef struct {
	hts_pos_t start, end;
} region1_t;

typedef struct bcf_sr_region_t {
	region1_t *regs;
	int nregs, mregs, creg;
} region_t;

const char *usage(void) {
	return
			"\n"
			"About: Extract CpG and nonCpG sites.\n"
			"Usage: mextr [file] [regions]\n"
			"Options:\n"
			"   -o, --cpgfile           Output file for CpG sites (default, not output)\n"
			"   -n, --noncpgfile        Output file for nonCpG sites (default, not output)\n"
			"   -b, --bed-methyl        Output file base for bedMethly files. Not compatible with multi-sample files  (default, not output)\n"
			"   -t, --bed-track-line    Track line for for bedMethly files (default, info taken from input VCF file)\n"
			"   -S, --report-file       Output file for JSON report (default, not output)\n"
			"   -r, --regions           restrict to comma separated list of regions\n"
			"   -R, --regions-file      restrict to regions listed in file\n"
			"   -@, --threads           Extra threads\n"
			"   -H, --no_header         Do not print header line(s) in output file(s) (default, false)\n"
			"   -g, --common-gt         Recall genotypes assuming common genotypes across samples\n"
			"   -m, --mode              Output mode for CpG sites\n"
			"         combined          Generate one line per CpG with combined estimates (default)\n"
			"         strand-specific   Generate two lines per CpG with the strand specific estimates\n"
			"   -w, --bw-mode           Output mode for bigWig files\n"
			"         combined          Generate one bigWig file for both strands (default)\n"
			"         strand-specific   Generate one bigWig files for each strand\n"
			"   -s, --select            Select mode for sites/CpGs\n"
			"         hom               Select on homozygote sites/CpGs (default)\n"
			"         het               Select on heterozygote sites/CpGs\n"
			"   -B, --reference-bias    Reference bias for re-calling (default 2)\n"
			"   -M, --min-nc            Minimum number of non-converted bases for non CpG site (default 1)\n"
			"   -p, --prop              Minimum proportion of sites/CpGs that must pass (default 0.0)\n"
			"   -N, --number            Minimum number of sites/CpGs that must pass (default 1)\n"
			"   -I, --inform            Minimum number of informative reads for a  CpG/site to pass (default 1)\n"
			"   -T, --threshold         Phred scaled threshold probability of selecting sites/CpGs (default 20)\n"
			"   -c, --conversion        <float>,<float> set under and over conversion rates\n"
			"   -Q, --bq-threshold      Base qality threshold used for calling\n"
			"   -z, --bgzip             Compress output with bgzip\n"
			"   -D, --md5               Calculate md5 digest for all output files\n"
			"   -x, --tabix             Generate tabix (tbx) indices for cpg and noncpg files\n"
			"\n";
}

static struct option loptions[] = {
		{"cpgfile",required_argument,0,'o'},
		{"noncpgfile",required_argument,0,'n'},
		{"bed-methyl",required_argument,0,'b'},
		{"bed-track-line",required_argument,0,'t'},
		{"report-file",required_argument,0,'S'},
		{"regions",required_argument,0,'r'},
		{"regions-file",required_argument,0,'R'},
		{"threads",required_argument,0,'@'},
		{"no_header",no_argument,0,'H'},
		{"common-gt",no_argument,0,'g'},
		{"mode",required_argument,0,'m'},
		{"bw-mode",required_argument,0,'w'},
		{"select",required_argument,0,'s'},
		{"prop",required_argument,0,'p'},
		{"min-nc",required_argument,0,'M'},
		{"reference-bias",required_argument,0,'B'},
		{"number",required_argument,0,'N'},
		{"inform",required_argument,0,'I'},
		{"threshold",required_argument,0,'T'},
		{"conversion",required_argument,0,'c'},
		{"bq-conversion",required_argument,0,'Q'},
		{"bgzip",no_argument,0,'z'},
		{"md5",no_argument,0,'D'},
		{"tabix",no_argument,0,'x'},
		{"help",no_argument,0,'h'},
		{0,0,0,0}
};

// Try to parse the paramaters used for bs_call from the headers
static void check_hdr_params(args_t *a) {
	char *par[] = {"under_conversion", "over_conversion", "mapq_thresh", "bq_thresh", NULL};
	bcf_hdr_t *h = a->hdr;
	for(int i = 0; i < h->nhrec; i++) {
		bcf_hrec_t *hr = h->hrec[i];
		if(hr->type == BCF_HL_GEN) {
			if(!strcmp(hr->key, "source") && !strncmp(hr->value, "bs_call", 7)) {
				char *p = strchr(hr->value, ',');
				while(p != NULL) {
					p++;
					int ix;
					for(ix = 0; par[ix] != NULL; ix++) if(!strncmp(p, par[ix], strlen(par[ix]))) break;
					if(par[ix] != NULL) {
						char *p1 = strchr(p, '=');
						if(p1) {
							switch(ix) {
							case 0:
								a->under_conv = strtod(p1 + 1, &p);
								break;
							case 1:
								a->over_conv = strtod(p1 + 1, &p);
								break;
							case 2:
								a->mq_thresh = (int)strtol(p1 + 1, &p, 10);
								break;
							case 3:
								a->bq_thresh = (int)strtol(p1 + 1, &p, 10);
								break;
							}
						}
					}
					p = strchr(p, ',');
				}
			}
		}
	}
}

static int cmp_reg(const void *s1, const void *s2, void *a) {
	const int *i1 = s1;
	const int *i2 = s2;
	const bcf_sr_regions_t *reg = a;
	return strcmp(reg->seq_names[*i1], reg->seq_names[*i2]);
}

struct tctg {
	char *name;
	uint64_t len;
};

static int cmp_tctg(const void *s1, const void *s2) {
	const struct tctg *ctg1, *ctg2;
	ctg1 = s1;
	ctg2 = s2;
	return strcmp(ctg1->name, ctg2->name);
}

void handle_command_line(int argc, char *argv[], args_t * const args) {
	int c;
	bool regions_file = false;
	char *regions_list = NULL;
	while ((c = getopt_long(argc, argv, "Q:Dxo:c:b:n:r:s:w:@:m:R:M:I:S:p:B:N:T:t:gzHah?",loptions,NULL)) >= 0) {
		switch (c) {
		case 'o':
			args->cpgfilename = optarg;
			break;
		case 'n':
			args->noncpgfilename = optarg;
			args->output_noncpg = true;
			break;
		case 'D':
			args->calc_md5 = true;
			break;
		case 'x':
			args->tabix = true;
			break;
		case 'S':
			args->reportfilename = optarg;
			break;
		case 'B':
			args->ref_bias = atof(optarg);
			break;
		case 'R':
			regions_file = true;
			// fall through
		case 'r':
			regions_list = optarg;
			break;
		case 'H':
			args->header = false;
			break;
		case 'g':
			args->common_gt = true;
			break;
		case 'w':
			if(!strcasecmp(optarg, "combined")) args->strand_specific = false;
			else if(!strcasecmp(optarg, "strand-specific")) args->strand_specific = true;
			break;
		case 's':
			if(!strcasecmp(optarg, "hom")) args->sel_mode = SELECT_HOM;
			else if(!strcasecmp(optarg, "het")) args->sel_mode = SELECT_HET;
			else error ("s (select) option can be either 'hom' or 'het'\n");
			break;
		case 'm':
			if(!strcasecmp(optarg, "combined")) args->mode = CPGMODE_COMBINED;
			else if(!strcasecmp(optarg, "strand-specific")) args->mode = CPGMODE_SEPARATE;
			else error ("m (mode) option can be either 'combined' or 'strand-specific'\n");
			break;
		case 'c':
			if (sscanf(optarg, "%lf,%lf", &args->under_conv, &args->over_conv) != 2)
				error("c (conversion) option expects two comma separated arguments)\n");
			break;
		case 'p':
			args->min_prop = atof(optarg);
			if(args->min_prop < 0.0) args->min_prop = 0.0;
			else if(args->min_prop > 1.0) args->min_prop = 1.0;
			break;
		case 'T':
			args->sel_thresh = atoi(optarg);
			if(args->sel_thresh < 0) args->sel_thresh = 0;
			else if(args->sel_thresh > 255) args->sel_thresh = 255;
			break;
		case 'N':
			args->min_num = atoi(optarg);
			if(args->min_num < 1) args->min_num = 1;
			break;
		case 'b':
			args->bedmethyl = optarg;
			break;
		case '@':
			args->threads = atoi(optarg);
			if(args->threads < 0) args->threads = 0;
			break;
		case 't':
			args->bedmethyl_track_line = optarg;
			break;
		case 'I':
			args->min_inform = atoi(optarg);
			if(args->min_inform < 0) args->min_inform = 0;
			break;
		case 'M':
			args->min_nc = atoi(optarg);
			if(args->min_nc < 0) args->min_nc = 0;
			break;
		case 'Q':
			args->bq_thresh = atoi(optarg);
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
	args->compress_threads = args->threads + 1;
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

	if(!bcf_sr_add_reader(args->sr, fname))
		error("failed to read from %s: %s\n", fname, bcf_sr_strerror(args->sr->errnum));

	args->hdr = args->sr->readers[0].header;
	if(args->sr->regions) {
		bcf_sr_regions_t *reg = args->sr->regions;
		// Sort regions by chromosome (required for bigWig, bigBed generation)
		int *ix = malloc(sizeof(int) * reg->nseqs);
		int nr = 0;
		for(int i = 0; i < reg->nseqs; i++) {
			if(bcf_hdr_name2id(args->hdr, reg->seq_names[i]) < 0) {
				fprintf(stderr,"Warning: requested contig %s not present in input file\n", reg->seq_names[i]);
				continue;
			}
			ix[nr++] = i;
		}
		if(!nr) error("None of the requested contigs are present in the input file\n");
		qsort_r(ix, nr, sizeof(int), cmp_reg, reg);
		char **tseq = malloc(sizeof(char *) * nr);
		region_t *treg = malloc(sizeof(region_t) * nr);
		for(int i = 0; i < nr; i++) {
			const int j = ix[i];
			tseq[i] = reg->seq_names[j];
			memcpy(treg + i, reg->regs + j, sizeof(region_t));
			khash_str2int_set(reg->seq_hash, tseq[i], i);
		}
		free(reg->seq_names);
		free(reg->regs);
		reg->seq_names = tseq;
		reg->regs = treg;
		reg->nseqs = nr;
		free(ix);
	}

	// If no regions have been specified.
	// extract all contigs from VCF/BCF header and add as regions after sorting
	if(!args->sr->regions) {
		int nctgs = args->hdr->n[BCF_DT_CTG];
		if(nctgs > 0) {
			// We can't add regions with an open reader so we remove (close) it, add the regions
			// and then re-open.  Removing the reader will result in the header info being
			// destroyed so we have to make a copy of the contig information beforehand
			struct tctg * const ctgs = malloc(sizeof(struct tctg) * nctgs);
			for(int i = 0; i < nctgs; i++) {
				const bcf_idpair_t * const idp = args->hdr->id[BCF_DT_CTG] + i;
				ctgs[i].name = strdup(idp->key);
				ctgs[i].len = idp->val->info[0];
			}
			qsort(ctgs, nctgs, sizeof(struct tctg), cmp_tctg);
			bcf_sr_destroy(args->sr);
			args->sr = bcf_sr_init();
			bcf_sr_set_threads(args->sr, args->threads);
			bcf_sr_regions_t * const reg = calloc(1, sizeof(bcf_sr_regions_t));
			reg->start = reg->end = -1;
			reg->prev_start = reg->prev_end = reg->prev_seq = -1;
			reg->seq_hash = khash_str2int_init();
			reg->seq_names = calloc(nctgs, sizeof(char *));
			reg->regs = calloc(nctgs, sizeof(region_t));
			for(int i = 0; i < nctgs; i++) {
				reg->nseqs++;
				reg->seq_names[i] = ctgs[i].name;
				khash_str2int_set(reg->seq_hash, reg->seq_names[i], i);
				reg->regs[i].creg = -1;
				reg->regs[i].nregs = reg->regs[i].mregs = 1;
				reg->regs[i].regs = malloc(sizeof(region1_t));
				reg->regs[i].regs->start = 0;
				reg->regs[i].regs->end = ctgs[i].len - 1;
			}
			args->sr->regions = reg;
			args->sr->explicit_regs = 1;
			args->sr->require_index = 1;
			free(ctgs);
			if(!bcf_sr_add_reader(args->sr, fname))
				error("failed to read from %s: %s\n", fname, bcf_sr_strerror(args->sr->errnum));
			args->hdr = args->sr->readers[0].header;
		}
	}
	bcf_sr_regions_t * const reg = args->sr->regions;
	int nctgs = reg->nseqs;
	args->cumul_len = malloc(sizeof(uint64_t) * nctgs);
	args->cumul_len[0] = reg->regs[0].regs->end + 1;
	for(int i = 1; i < nctgs; i++) args->cumul_len[i] = args->cumul_len[i - 1] + reg->regs[i].regs->end + 1;
	int nctgs_vcf = args->hdr->n[BCF_DT_CTG];
	args->id_trans = malloc(sizeof(int) * nctgs_vcf);
	for(int i = 0; i < nctgs_vcf; i++) {
		const bcf_idpair_t * const idp = args->hdr->id[BCF_DT_CTG] + i;
		int ret = khash_str2int_get(args->sr->regions->seq_hash, idp->key, args->id_trans + i);
		if(ret < 0) args->id_trans[i] = -1;
	}
	int ns = bcf_hdr_nsamples(args->hdr);
	assert(ns > 0);
	for(int i = 0; i < REC_BUF_SIZE; i++) args->rec_buf.buf[i] = rec_init(ns);
	if((args->bedmethyl) && ns > 1) error("bedMethyl output not compatible with multi-sample files\n");
	check_hdr_params(args);
}
