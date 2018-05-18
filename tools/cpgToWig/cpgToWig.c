#include "config.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <signal.h>
#include <locale.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>

#include "libhdr.h"
#include "utils.h"
#include "lk_malloc.h"
#include "loki_compress.h"
#include "string_utils.h"
#include "lkgetopt.h"

#define DEFAULT_PHRED (20)
#define DEFAULT_MIN_INF (5)

static void usage(FILE *f) {
  fputs("usage:\n cpgToWig <input file> \n", f);
  fputs("  -o|--out_prefix     Prefix for the output file names\n",f);
  fprintf(f, "  -q|--quality        PHRED threshold for genotype calling       (default='%d')\n", DEFAULT_PHRED);
  fprintf(f, "  -m|--min_inform     Minimum informative reads                  (default='%d')\n", DEFAULT_MIN_INF);
  fputs("  -h|help|usage       Print this help \n\n", f);
}

typedef struct {
  char *out_prefix;
  int threshold;
  int min_inf;
} filter_params;

static int process_file(FILE *fp, filter_params *params) {
  FILE *fp_call = NULL, *fp_cov = NULL;
  assert(params->out_prefix != NULL);
  size_t l = strlen(params->out_prefix);
  char *outfile = lk_malloc(l + 10);
  memcpy(outfile, params->out_prefix, l);
  strcpy(outfile + l, "_call.wig");
  fp_call = fopen(outfile, "w");
  if(fp_call == NULL) return errno;
  strcpy(outfile + l, "_cov.wig");
  fp_cov = fopen(outfile, "w");
  if(fp_cov == NULL) return errno;
  free(outfile);

  tokens *tok = NULL;
  string *s = NULL;
  void *tbuf = NULL;
  char *prev_chrom = NULL;
  int err = 0;
  while (!err) {
    s = fget_string(fp, s, &tbuf);
    if (!s->len) break;
    tok = tokenize(get_cstring(s), '\t', tok);
    if(tok->n_tok == 13) {
      if(!strcmp(tok->toks[3], "CG") && tok->toks[5][0] != '-') {
	int qual = atoi(tok->toks[4]);
	int non_conv = atoi(tok->toks[7]);
	int conv = atoi(tok->toks[8]);
	int cov = non_conv + conv;
	if(qual >= params->threshold && cov >= params->min_inf) {
	  if(!prev_chrom || strcmp(prev_chrom, tok->toks[0])) {
	    fprintf(fp_call, "variableStep\tchrom=%s\tspan=2\n", tok->toks[0]);
	    fprintf(fp_cov, "variableStep\tchrom=%s\tspan=2\n", tok->toks[0]);
	    if(prev_chrom) free(prev_chrom);
	    size_t l = strlen(tok->toks[0]) + 1;
	    prev_chrom = lk_malloc(l);
	    memcpy(prev_chrom, tok->toks[0], l);
	  }
	  fprintf(fp_call, "%s\t%s\n", tok->toks[1], tok->toks[5]);
	  fprintf(fp_cov, "%s\t%d\n", tok->toks[1], cov);
	}
      }
    }
  }

  fclose(fp_cov);
  fclose(fp_call);
  if(s) free_string(s);
  if(tok) free_tokens(tok);
  if(tbuf) free_fget_buffer(&tbuf);

  return err;
}

int main(int argc, char *argv[]) {

  filter_params params = {
    .threshold = DEFAULT_PHRED,
    .min_inf = DEFAULT_MIN_INF,
    .out_prefix = NULL,
  };

  static struct option longopts[] = {
    {"output_prefix", required_argument, 0, 'o'},
    {"quality", required_argument, 0, 'q'},
    {"min_informative", required_argument, 0, 'm'},
    {"help", no_argument, 0, 'h'},
    {"usage", no_argument, 0, 'h'},
    {0, 0, 0, 0}};

  int err = 0;
  int c;
  while (!err && (c = getopt_long(argc, argv, "o:q:m:h?", longopts, 0)) != -1) {
    switch (c) {
    case 'o':
      params.out_prefix = optarg;
      break;
    case 'q':
      params.threshold = atoi(optarg);
      break;
    case 'm':
      params.min_inf= atoi(optarg);
      break;
    case 'h':
    case '?':
      err = 1;
      usage(stdout);
      break;
    }
  }
  if (err == 1) return 0;
  if(params.out_prefix == NULL) params.out_prefix = "cpgToWig";
  FILE *fp;
  if (argc > optind) {
    int j;
    fp = open_readfile_and_check(argv[optind], &j, init_compress());
  } else fp = stdin;
  if (fp) err = process_file(fp, &params);
  if (fp != stdin) fclose(fp);
  return err;
}
