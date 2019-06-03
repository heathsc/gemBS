#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <sys/wait.h>

#include "mextr.h"
#include "utils.h"
#include "compress.h"

static FILE *open_ofile(char *name, int compress, bool append) {
  FILE *fp = NULL;
  if(append) compress = 0;
  int comp_ix = COMPRESS_NONE;
  if(compress != 0) {
    for(comp_ix = 0; comp_ix < COMPRESS_NONE; comp_ix++) {
      if(compress & (1 << comp_ix)) break;
    }
  }
  struct compress *cdata = comp_ix < COMPRESS_NONE ? get_compress_data() : NULL;
  if(name != NULL) {
    if(comp_ix < COMPRESS_NONE) {
      char *tname = name;
      char *suffix = cdata->compress_suffix[comp_ix];
      // Check whether file name already has suffix
      char *p = strrchr(tname, '.');
      if(p == NULL || strcmp(p + 1, suffix)) {
	// No, so we will have to add it
	tname = malloc(strlen(name) + strlen(suffix) + 2);
	sprintf(tname, "%s.%s", name, suffix);
      }
      int i = child_open(WRITE, tname, cdata->comp_path[comp_ix][0]);
      fp = fdopen(i, "w");
      if(name != tname) free(tname);
    } else {
      fp = append ? fopen(name, "a") : fopen(name, "w");
    }
  } else {
    if(isatty(fileno(stdout))) comp_ix = COMPRESS_NONE;
    if(comp_ix < COMPRESS_NONE) {
      int i = child_open(WRITE, NULL, cdata->comp_path[comp_ix][0]);
      fp = fdopen(i, "w");
    } else fp = stdout;
  }
  return fp;
}

static void init_stats(args_t *a) {
  a->stats = calloc((size_t)1, sizeof(stats_t));
}

static void write_stats(args_t *a) {
  if(a->stats != NULL) {
    a->reportfile = a->reportfilename == NULL ? NULL : fopen(a->reportfilename, "w");
    if(a->reportfile != NULL) {
      FILE *fp = a->reportfile;
      stats_t *st = a->stats;
      fprintf(fp,"{\n\t\"TotalSites\": %" PRIu64 ",\n", st->n_sites);
      fprintf(fp,"\t\"SitesPassed\": %" PRIu64 "\n", st->n_sites_pass);
      fputs("}\n", fp);
      fclose(fp);
    }
  }
}

static gt_meth *sample_gt[2];
static cpg_prob *sample_cpg;
static double *sample_Q[3];

static void init_files(args_t *a) {
  if(a->cpgfilename == NULL) a->cpgfile = NULL;
  else if(a->cpgfilename[0] == '-' && a->cpgfilename[1] == 0) a->cpgfile = open_ofile(NULL, 0, false);
  else  a->cpgfile = open_ofile(a->cpgfilename, a->compress, a->append_mode);
  if(a->wigfilename == NULL) a->wigfile = NULL;
  else if(a->wigfilename[0] == '-' && a->wigfilename[1] == 0) a->wigfile = open_ofile(NULL, 0, false);
  else  a->wigfile = open_ofile(a->wigfilename, a->compress, a->append_mode);
  
  a->noncpgfile = a->noncpgfilename == NULL ? (a->output_noncpg ? a->cpgfile : NULL) : open_ofile(a->noncpgfilename, a->compress, a->append_mode);
  if(a->bedmethyl != NULL) {
    char *p = strrchr(a->bedmethyl, '.');
    if(p && !strcmp(".bed",p)) {
      *p = 0;
    }
    p = strrchr(a->bedmethyl, '_');
    if(p && !strcmp("_cpg",p)) {
      *p = 0;
    }
    size_t l = strlen(a->bedmethyl);
    p = malloc((l + 9) * 3);
    a->bedmethylnames[0] = p;
    a->bedmethylnames[1] = p + l + 9;
    a->bedmethylnames[2] = p + 2 * (l + 9);
    sprintf(a->bedmethylnames[BEDMETHYL_CPG], "%s_cpg.bed", a->bedmethyl);
    sprintf(a->bedmethylnames[BEDMETHYL_CHG], "%s_chg.bed", a->bedmethyl);
    sprintf(a->bedmethylnames[BEDMETHYL_CHH], "%s_chh.bed", a->bedmethyl);
    for(int i = 0; i < 3; i++) 
      a->bedmethylfiles[i] = open_ofile(a->bedmethylnames[i], a->compress, a->append_mode);
  }
}

static void print_file_header(FILE *fp, int ns, char **names) {
  if(fp != NULL) {
    fputs("Contig\tPos0\tPos1\tRef", fp);
    for(int i = 0; i < ns; i++) {
      char *name = names[i];
      fprintf(fp, "\t%s:Call\t%s:Flags\t%s:Meth\t%s:non_conv\t%s:conv\t%s:support_call\t%s:total", name, name, name, name, name, name, name);
    }
		fputc('\n', fp);
  }
}

char *copy_and_strip_quotes(char *s) {
  if(!s) return s;
  size_t l = strlen(s);
  if(l > 1) {
    if((s[0] == '\"' && s[l-1] =='\"') || (s[0] == '\'' && s[l-1] =='\'')) {
      (s++)[--l] = 0;
    }
  }
  char *s1 = malloc(l + 1);
  if(s1 != NULL) memcpy(s1, s, l + 1);
  return s1;
}

static void print_bedmethyl_headers(args_t *args) {
  if(args->bedmethyl_track_line == NULL) {
    char *sample_name = NULL;
    char *sample_desc = NULL;
    char *sample_bc = NULL;
    // Try and get sample info from VCF file headers
    bcf_hdr_t *h = args->hdr;
    for(int i = 0; i < h->nhrec; i++) {
      bcf_hrec_t *hr = h->hrec[i];
      if(hr->type == BCF_HL_STR) {
	if(!strcmp(hr->key, "bs_call_sample_info")) {
	  int ix = bcf_hrec_find_key(hr, "ID");
	  if(ix >= 0) {
	    sample_bc = copy_and_strip_quotes(hr->vals[ix]);
	    ix = bcf_hrec_find_key(hr, "SM");
	    if(ix >= 0) sample_name = copy_and_strip_quotes(hr->vals[ix]);
	    ix = bcf_hrec_find_key(hr, "DS");
	    if(ix >= 0) sample_desc = copy_and_strip_quotes(hr->vals[ix]);
	  }
	}
      }
    }
    if(sample_name == NULL) sample_name = strdup(h->samples[0]);
    if(sample_desc == NULL) sample_desc = strdup(sample_name);
    for(bedmethyl_type t = BEDMETHYL_CPG; t <= BEDMETHYL_CHH; t++) {
      FILE *fp = args->bedmethylfiles[t];
      if(fp != NULL) {
	fprintf(fp, "track name=\"%s\" description=\"%s\" visibility=2 itemRgb=\"On\"\n", sample_desc, sample_name);
      }
    }
    if(sample_bc) free(sample_bc);
    if(sample_name) free(sample_name);
    args->bedmethyl_desc = sample_desc;
  } else {
    for(bedmethyl_type t = BEDMETHYL_CPG; t <= BEDMETHYL_CHH; t++) {
      char *line = args->bedmethyl_track_line;
      size_t l = strlen(line);
      if(l > 1 && line[l - 1] == '\n') line[--l] = 0;
      if(!strncmp(line, "track ", 6)) line += 6;
      FILE *fp = args->bedmethylfiles[t];
      if(fp != NULL) fprintf(fp, "track %s\n", line);
    }
  }
}

static void print_headers(args_t *args) {
  int ns = bcf_hdr_nsamples(args->hdr);
  if(args->cpgfile) print_file_header(args->cpgfile, ns, args->hdr->samples);
  if(args->output_noncpg && args->noncpgfile != args->cpgfile)
    print_file_header(args->noncpgfile, ns, args->hdr->samples);
  print_bedmethyl_headers(args);
}

static void close_files(args_t *a) {
  if(a->cpgfile != NULL && a->cpgfile != stdout) fclose(a->cpgfile);
  if(a->noncpgfile != NULL) fclose(a->noncpgfile);
  if(a->wigfile != NULL && a->wigfile != stdout) fclose(a->wigfile);
  for(int i = 0; i < 3; i++)
    if(a->bedmethylfiles[i] != NULL) fclose(a->bedmethylfiles[i]);
  while(waitpid(-1, NULL, 0) > 0);
}

static args_t args = {
  .hdr = NULL,
  .cpgfile = NULL,
  .noncpgfile = NULL,
  .wigfile = NULL,
  .reportfile = NULL,
  .cpgfilename = NULL,
  .wigfilename = NULL,
  .bedmethylfiles = {NULL, NULL, NULL},
  .noncpgfilename = NULL,
  .reportfilename = NULL,
  .bedmethyl = NULL,
  .bedmethylnames = {NULL, NULL, NULL},
  .bedmethyl_track_line = NULL,
  .bedmethyl_desc = ".",
  .stats = NULL,
  .min_prop = 0.0,
  .min_num = 1,
  .min_inform = 0,
  .min_nc = 1,
  .ref_bias = DEFAULT_REF_BIAS,
  .under_conv = DEFAULT_UNDER_CONV,
  .over_conv = DEFAULT_OVER_CONV,
  .bq_thresh = DEFAULT_BQ_THRESH,
  .mq_thresh = DEFAULT_MAPQ_THRESH,
  .mode = CPGMODE_COMBINED,
  .sel_mode = SELECT_HOM,
  .sel_thresh = DEFAULT_SELECT_THRESH,
  .compress = 0,
  .common_gt = false,
  .output_noncpg = false,
  .header = true,
  .append_mode = false
};

const char *about(void)
{
  return "Extract CpG and nonCpG sites.\n";
}

const char *usage(void)
{
  return
    "\n"
    "About: Extract CpG and nonCpG sites.\n"
    "Usage: bcftools +mextr [General Options] -- [Plugin Options]\n"
    "Options:\n"
    "   run \"bcftools plugin\" for a list of common options\n"
    "\n"
    "Plugin options:\n"
    "   -o, --cpgfile           Output file for CpG sites (default = stdout)\n"
    "   -n, --noncpgfile        Output file for nonCpG sites (default, not output)\n"
    "   -b. --bed-methyl        Output file base for bedMethly files. Not compatible with multi-sample files  (default, not output)\n" 
    "   -w. --wigfile           Output file for wig file (methylation)\n" 
    "   -t. --bed-track-line    Track line for for bedMethly files (default, info taken from input VCF file)\n" 
    "   -r, --report-file       Output file for JSON report (default, not output)\n"
    "   -H, --no_header         Do not print header line(s) in output file(s) (default, false)\n"
    "   -g, --common-gt         Recall genotypes assuming common genotypes across samples\n"
    "   -m, --mode              Output mode for CpG sites\n"
    "         combined          Generate one line per CpG with combined estimates (default)\n" 
    "         strand-specific   Generate two lines per CpG with the strand specific estimates\n"
    "   -s, --select            Select mode for sites/CpGs\n"
    "         hom               Select on homozygote sites/CpGs (default)\n"
    "         het               Select on heterozygote sites/CpGs\n"
    "   -R, --reference-bias    Reference bias for re-calling (default 2)\n"
    "   -M, --min-nc            Minimum number of non-converted bases for non CpG site (default 1)\n"
    "   -p, --prop              Minimum proportion of sites/CpGs that must pass (default 0.0)\n"
    "   -N, --number            Minimum number of sites/CpGs that must pass (default 1)\n"
    "   -I, --inform            Minimum number of informative reads for a  CpG/site to pass (default 1)\n"
    "   -T, --threshold         Phred scaled threshold probability of selecting sites/CpGs (default 20)\n"
    "   -c, --conversion        <float>,<float> set under and over conversion rates\n"
    "   -Q, --bq-threshold      Base qality threshold used for calling\n"
    "   -z, --gzip              Compress output with gzip (bgzip if available)\n"
    "   -j, --bzip2             Compress output with bzip2\n"
    "   -x, --xz                Compress output with xz\n"
    "   -a, --append            Append to output files rather than create new ones.  Not compatible wih output compression\n"
    "\n"
    "Example:\n"
    "   bcftools +mextr in.vcf -- -o out_cpg.txt -n out_noncpg.txt -z\n"
    "\n";
}

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

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out __unused__)
{
  args.hdr  = in;
  check_hdr_params(&args);
  static struct option loptions[] = {
    {"cpgfile",required_argument,0,'o'},
    {"wigfile",required_argument,0,'w'},
    {"noncpgfile",required_argument,0,'n'},
    {"bed-methyl",required_argument,0,'b'},
    {"bed-track-line",required_argument,0,'t'},
    {"report-file",required_argument,0,'r'},
    {"no_header",no_argument,0,'H'},
    {"common-gt",no_argument,0,'g'},
    {"mode",required_argument,0,'m'},
    {"select",required_argument,0,'s'},
    {"prop",required_argument,0,'p'},
    {"min-nc",required_argument,0,'M'},
    {"reference-bias",required_argument,0,'R'},
    {"number",required_argument,0,'N'},
    {"inform",required_argument,0,'I'},
    {"threshold",required_argument,0,'T'},
    {"conversion",required_argument,0,'c'},
    {"bq-conversion",required_argument,0,'Q'},
    {"gzip",no_argument,0,'z'},
    {"bzip2",no_argument,0,'j'},
    {"xz",no_argument,0,'x'},
    {"append",no_argument,0,'a'},
    {0,0,0,0}
  };
  int c;
  bool mult_comp = false;
  while ((c = getopt_long(argc, argv, "?Qh:o:c:b:n:r:m:R:M:I:s:p:N:T:t:w:gzHjxa",loptions,NULL)) >= 0) {
    switch (c) {
    case 'a':
      args.append_mode = true;
      break;
    case 'o':
      args.cpgfilename = optarg;
      break;
    case 'w':
      args.wigfilename = optarg;
      break;
    case 'n':
      args.noncpgfilename = optarg;
      args.output_noncpg = true;
      break;
    case 'r':
      args.reportfilename = optarg;
      break;
    case 'R':
      args.ref_bias = atof(optarg);
      break;
    case 'H':
      args.header = false;
      break;
    case 'g':
      args.common_gt = true;
      break;
    case 's':
      if(!strcasecmp(optarg, "hom")) args.sel_mode = SELECT_HOM;
      else if(!strcasecmp(optarg, "het")) args.sel_mode = SELECT_HET;
      else error ("s (select) option can be either 'hom' or 'het'\n");
      break;
    case 'm':
      if(!strcasecmp(optarg, "combined")) args.mode = CPGMODE_COMBINED;
      else if(!strcasecmp(optarg, "strand-specific")) args.mode = CPGMODE_SEPARATE;
      else error ("m (mode) option can be either 'combined' or 'strand-specific'\n");
      break;
    case 'c':
      if (sscanf(optarg, "%lf,%lf", &args.under_conv, &args.over_conv) != 2) 
	error("c (conversion) option expects two comma separated arguments)\n");
      break;
    case 'p':
      args.min_prop = atof(optarg);
      if(args.min_prop < 0.0) args.min_prop = 0.0;
      else if(args.min_prop > 1.0) args.min_prop = 1.0;
      break;
    case 'T':
      args.sel_thresh = atoi(optarg);
      if(args.sel_thresh < 0) args.sel_thresh = 0;
      else if(args.sel_thresh > 255) args.sel_thresh = 255;
      break;
    case 'N':
      args.min_num = atoi(optarg);
      if(args.min_num < 1) args.min_num = 1;
      break;
    case 'b':
      args.bedmethyl = optarg;
      break;
    case 't':
      args.bedmethyl_track_line = optarg;
      break;
    case 'I':
      args.min_inform = atoi(optarg);
      if(args.min_inform < 0) args.min_inform = 0;
      break;
    case 'M':
      args.min_nc = atoi(optarg);
      if(args.min_nc < 0) args.min_nc = 0;
      break;
    case 'Q':
      args.bq_thresh = atoi(optarg);
      break;
    case 'z':
      if(args.compress) mult_comp = true;
      else args.compress |= COMP_GZIP;
      break;
    case 'j':
      if(args.compress) mult_comp = true;
      else args.compress |= COMP_BZIP2;
      break;
    case 'x':
      if(args.compress) mult_comp = true;
      else args.compress |= COMP_XZ;
      break;
    case 'h':
    case '?':
    default: error(usage()); break;
    }
  }
  if(mult_comp) error("Can not combine multiple compression options\n");
  if(args.append_mode && args.compress) error("Output compression not compatible with append mode\n");
  int ns = bcf_hdr_nsamples(args.hdr);
  assert(ns > 0);
  if((args.bedmethyl || args.wigfile) && ns > 1) error("bedMethyl and wig output not compatible with multi-sample files\n");
  if (optind != argc) error(usage());
  init_files(&args);
  if(!args.append_mode && (args.header || args.bedmethyl)) print_headers(&args);
  if(args.reportfilename != NULL) init_stats(&args);
  sample_gt[0] = malloc(sizeof(gt_meth) * ns * 2);
  sample_gt[1] = sample_gt[0] + ns;
  sample_Q[0] = malloc(sizeof(double) * (ns + 1) * 2 + ns);
  sample_Q[1] = sample_Q[0] + ns + 1;
  sample_Q[2] = sample_Q[1] + ns + 1;
  if(args.mode == CPGMODE_COMBINED) sample_cpg = malloc(sizeof(cpg_prob) * ns);
  fill_base_prob_table();
  return 1;
}

static fmt_field_t tags[] = {
			     { "FT", BCF_HT_STR, {{NULL, 0, 0}, {NULL, 0, 0}}},
			     { "MC8", BCF_HT_INT, {{NULL, 0, 0}, {NULL, 0, 0}}},
			     { "AMQ", BCF_HT_INT, {{NULL, 0, 0}, {NULL, 0, 0}}},
			     { "CX", BCF_HT_STR, {{NULL, 0, 0}, {NULL, 0, 0}}},
			     { "AQ", BCF_HT_INT, {{NULL, 0, 0}, {NULL, 0, 0}}},
			     { "MQ", BCF_HT_INT, {{NULL, 0, 0}, {NULL, 0, 0}}},
			     { NULL, 0, {{NULL, 0, 0}, {NULL, 0, 0}}},
};

bcf1_t *process(bcf1_t *rec)
{
  static int idx;
  static int32_t curr_rid = -1, prev_pos;
  static bool valid[2] = {false, false};
  static bcf1_t prev_rec;
  
  int ns = bcf_hdr_nsamples(args.hdr);
  stats_t *st = args.stats;
  if(st != NULL) st->n_sites++;
  bcf_unpack(rec, BCF_UN_FLT);
  int n_all = rec->n_allele;
  bool cg = false;
  for(int i = 0; i < n_all; i++) {
    char c = rec->d.allele[i][0];
    if((c == 'C' || c == 'G') && rec->d.allele[i][1] == 0) {
      cg = true;
      break;
    }
  }
  if(cg || args.output_noncpg) { // Site with potentially Cs or Gs (or we are outputting non_cpgs)
    bcf_unpack(rec, BCF_UN_ALL);
    // Get format tags
    for(int ix = 0; tags[ix].tag != NULL; ix++) {
      fmt_store_t *s = tags[ix].st + idx;
      s->ne = bcf_get_format_values(args.hdr, rec, tags[ix].tag, &s->dat_p, &s->dat_n, tags[ix].type);
    }
    if(tags[FMT_CX].st[idx].ne > 0 && tags[FMT_MC8].st[idx].ne == ns * 8) {
      // Get sample base counts and genotype probs.
      int32_t *mc8_p = tags[FMT_MC8].st[idx].dat_p;
      int32_t *amq_p = tags[FMT_AMQ].st[idx].dat_p;
      int n_amq = tags[FMT_AMQ].st[idx].ne / ns;
      int32_t *aq_p = tags[FMT_AQ].st[idx].ne == ns ? tags[FMT_AQ].st[idx].dat_p : NULL; 
      int32_t *mq_p = tags[FMT_MQ].st[idx].ne == ns ? tags[FMT_MQ].st[idx].dat_p : NULL;
      double ms_mq = 0.0;
      int32_t tot_n = 0;
      for(int i = 0; i < ns; i++) {
	int32_t *ct = sample_gt[idx][i].counts;
	int32_t *amq = sample_gt[idx][i].aqual;
	memset(ct, 0, sizeof(int32_t) * 8);
	memset(amq, 0, sizeof(int32_t) * 8);
	int32_t x = mc8_p[i * 8];
	int k = 0;
	if(x != bcf_int32_missing) {
	  int k1 = 0;
	  for(int j = 0; j < 8; j++) {
	    x = mc8_p[i * 8 + j];
	    ct[j] += x;
	    k += x;
	    if(x > 0 && amq_p != NULL && k1 < n_amq) {
	      int q = amq_p[i * n_amq + k1++];
	      if(q >= 0) {
		if(q > MAX_QUAL) q = MAX_QUAL;
		amq[j] = q;
	      }
	    }
	  }
	  if(amq_p == NULL) {
	    int q = aq_p == NULL ? args.bq_thresh : aq_p[i];
	    if(q > MAX_QUAL) q = MAX_QUAL;
	    for(int j = 0; j < 8; j++) amq[j] = q;
	  }
	}
	if(k > 0) {	  
	  if(mq_p != NULL) {
	    int m = mq_p[i];
	    ms_mq += (double)k * (double)(m * m);
	  }
	  tot_n += k;
	  calc_gt_prob(sample_gt[idx] + i, &args, rec->d.allele[0][0]);
	  sample_gt[idx][i].skip = false;
	} else sample_gt[idx][i].skip = true;
      }
      // If we force a common genotype, calculate prob. distribution for common genotype
      if(args.common_gt) {
	double gt[10];
	for(int k = 0; k < 10; k++) gt[k] = 0.0;
	for(int i = 0; i < ns; i++) {
	  if(!sample_gt[idx][i].skip) {
	    for(int k = 0; k < 10; k++) gt[k] += sample_gt[idx][i].gt_prob[k];
	  }
	}
	double max = gt[0];
	int max_gt = 0;
	for(int k = 1; k < 10; k++) {
	  if(gt[k] > max) {
	    max = gt[k];
	    max_gt = k;
	  }
	}
	double sum = 0.0;
	for(int k = 0; k < 10; k++) sum += exp(gt[k] - max);
	sum = log(sum);
	for(int k = 0; k < 10; k++) gt[k] -= (max + sum);
	for(int i = 0; i < ns; i++) {
	  if(!sample_gt[idx][i].skip) {
	    for(int k = 0; k < 10; k++) sample_gt[idx][i].gt_prob[k] = gt[k];
	    sample_gt[idx][i].max_gt = max_gt;
	    sample_gt[idx][i].sum = max + sum;
	  }
	}
      }
      valid[idx] = true;
      // Here is the logic for deciding what we print
      if(rec->rid != curr_rid) curr_rid = rec->rid;
      else if(rec->pos - prev_pos == 1 && valid[idx ^ 1]) output_cpg(&args, &prev_rec, tags, sample_gt, idx ^ 1, sample_cpg, sample_Q);
      if(args.bedmethyl || args.wigfile) {
	output_bedmethyl(&args, rec, tags, sample_gt, idx);
      }
      idx ^= 1;
      prev_pos = rec->pos;
      memcpy(&prev_rec, rec, sizeof(bcf1_t));
      valid[idx] = false;
      if(st != NULL) st->n_sites_pass++;
    }
  }
  return NULL;
}

void destroy(void)
{
  write_stats(&args);
  close_files(&args);
}
