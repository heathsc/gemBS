#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <math.h>
#include <stdbool.h>

#include "snpxtr.h"
#include "utils.h"
#include "compress.h"

static args_t args = {
  .hdr = NULL,
  .output_filename = NULL,
  .snp_filename = NULL,
  .output_file = NULL,
  .snp_hash = NULL,
  .gt = NULL,
  .compress = 0,
  .pass_index = -1,
  .gt_index = -1,
  .dbSNP = NULL,
  .dbSNP_header = NULL,
  .dbSNP_prefix = NULL,
  .dbSNP_name = NULL,
  .n_dbSNP_prefixes = 0
};

const char *about(void)
{
  return "Extract SNPs\n";
}

const char *usage(void)
{
  return
    "\n"
    "About: Extract SNPs from VCF file.\n"
    "Usage: bcftools +snpxtr [General Options] -- [Plugin Options]\n"
    "Options:\n"
    "   run \"bcftools plugin\" for a list of common options\n"
    "\n"
    "Plugin options:\n"
    "   -o, --output            Output file (default = stdout)\n"
    "   -s, --snps              File containing list of SNP IDs to be selected (default selected all sites with PASS))\n"
    "   -D, --dbsnp             (dbSNP processed file)\n"
    "   -z, --gzip              Compress output with gzip (bgzip if available)\n"
    "   -j, --bzip2             Compress output with bzip2\n"
    "   -x, --xz                Compress output with xz\n"
    "\n"
    "Example:\n"
    "   bcftools +snpxtr in.vcf -- -s snp_list.txt -o out_snps.txt -z\n"
    "\n";
}

static FILE *open_ofile(char *name, int compress) {
  FILE *fp = NULL;
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
	asprintf(&tname, "%s.%s", name, suffix);
      }
      int i = child_open(WRITE, tname, cdata->comp_path[comp_ix][0]);
      fp = fdopen(i, "w");
      if(name != tname) free(tname);
    } else fp = fopen(name, "w");
  } else {
    if(isatty(fileno(stdout))) comp_ix = COMPRESS_NONE;
    if(comp_ix < COMPRESS_NONE) {
      int i = child_open(WRITE, NULL, cdata->comp_path[comp_ix][0]);
      fp = fdopen(i, "w");
    } else fp = stdout;
  }
  return fp;
}

static void init_files(args_t *a) {
  a->output_file = open_ofile(a->output_filename, a->compress);
}

static void close_files(args_t *a) {
  if(a->output_filename != NULL && a->output_file != stdout) fclose(a->output_file);
}

void read_snp_file(void) {
  bool filter;
  FILE *fp = open_readfile_and_check(args.snp_filename, &filter);
  char *buf = NULL;
  size_t buf_size = 0;
  tokens *tok = NULL;
  int nsnps = 0;
  fprintf(stderr,"Reading SNP list from %s\n", args.snp_filename);
  for(;;) {
    ssize_t l = getline(&buf, &buf_size, fp);
    if(l < 0) break;
    tok = tokenize(buf, '\t', tok);
    if(tok->n_tok >= 1) {
      char *p = tok->toks[0];
      char *id = malloc(strlen(p) + 3);
      sprintf(id, "rs%s", p);
      snp *s;
      HASH_FIND_STR(args.snp_hash, id, s);
      if(s == NULL) {
	nsnps++;
	s = malloc(sizeof(snp));
	s->name = id;
	HASH_ADD_KEYPTR(hh, args.snp_hash, s->name, strlen(s->name), s);
      }      
    }
  }
  fclose(fp);
  if(buf) free(buf);
  if(filter) {
    int i;
    while(waitpid(-1, &i, WNOHANG) > 0);
  }
  fprintf(stderr,"%d SNPs read in\n", nsnps);
}

static void store_dbsnp_entries(dbsnp_bin *bin, int n_entries, int name_buf_sz, uint16_t *entries, uint8_t *name_buf) {
  bin->entries = malloc(sizeof(uint16_t) * n_entries);
  bin->name_buf = malloc((size_t)name_buf_sz);
  bin->n_entries = n_entries;
  uint64_t msk = (uint64_t)0;
  for(int i = 0; i < n_entries; i++) {
    bin->entries[i] = entries[i];
    msk |= ((uint64_t)1 << (entries[i] & 63));
  }
  bin->mask = msk;
  memcpy(bin->name_buf, name_buf, name_buf_sz);
}

static uint8_t db_tab[] = {
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
  0xff, 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x10, 0x11, 0x12, 0x13, 0x14, 
  0x15, 0x16, 0x17, 0x18, 0x19, 0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x30,
  0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46,
  0x47, 0x48, 0x49, 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x60, 0x61, 0x62,
  0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78,
  0x79, 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x90, 0x91, 0x92, 0x93, 0x94,
  0x95, 0x96, 0x97, 0x98, 0x99, 0x0f, 0x1f, 0x2f, 0x3f, 0x4f, 0x5f, 0x6f, 0x7f, 0x8f, 0x9f, 0xff,
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff
};
  
static void read_dbSNP_file(void) {
  bool ok = true;
  bool filter;
  FILE *fp = open_readfile_and_check(args.dbSNP_name, &filter);
  char *buf = NULL;
  size_t buf_size = 0;
  fprintf(stderr,"Loading dbSNP from %s\n", args.dbSNP_name);
  int cbin = 0, max_bin = 0, n_entries = 0, name_buf_ptr = 0, ix = -1;
  dbsnp_bin *bins = NULL;
  char *ctg_name = NULL;
  dbsnp_ctg *ctg = NULL;
  uint16_t *entries = malloc(sizeof(uint16_t) * 64);
  uint8_t *name_buf = malloc(sizeof(uint8_t) * 256 * 64);
  int n_snps = 0, n_bins = 0, n_ctgs = 0;
  ssize_t l = getline(&buf, &buf_size, fp);
  if(l > 0 && buf[l - 1] == '\n') buf[--l] = 0;
  if(l > 7 && !strncmp(buf, "track ", 6)) {
    args.dbSNP_header = malloc((size_t)(l - 5));
    memcpy(args.dbSNP_header, buf + 6, l - 6);
    args.dbSNP_header[l - 6] = 0;
    int n_p_store = 8;
    char **p_store = malloc(sizeof(void *) * n_p_store);
    while(true) {
      l = getline(&buf, &buf_size, fp);
      if(l > 0 && buf[l - 1] == '\n') buf[--l] = 0;
      if(l < 1 || buf[0] != '+') break;
      if(l > 1) {
	if(args.n_dbSNP_prefixes == 0xffff) {
	  fprintf(stderr, "Error in dbSNP file: too many prefixes\n");
	  exit(-1);
	}
	if(args.n_dbSNP_prefixes == n_p_store) {
	  n_p_store *= 1.5;
	  p_store = realloc(p_store, sizeof(void *) * n_p_store);
	}
	char *tp = p_store[args.n_dbSNP_prefixes++] = malloc(l);
	memcpy(tp, buf + 1, l - 1);
	tp[l - 1] = 0;
      }
    }
    if(!args.n_dbSNP_prefixes) {
      fprintf(stderr, "Error in dbSNP file: no prefix information\n");
      ok = false;
    } else {
      args.dbSNP_prefix = malloc(sizeof(void *) * args.n_dbSNP_prefixes);
      memcpy(args.dbSNP_prefix, p_store, sizeof(void *) * args.n_dbSNP_prefixes);
      free(p_store);
    }
    while(l >= 0 && ok) {
      if(l > 0) {
	if(buf[0] == '>') {
	  if(n_entries) {
	    store_dbsnp_entries(bins, n_entries, name_buf_ptr, entries, name_buf);
	    n_bins++;
	    n_snps += n_entries;
	  }
	  if(cbin != max_bin) {
	    fprintf(stderr, "Error in dbSNP file - wrong number of bins (expected %d, saw %d\n", max_bin, cbin);
	    ok = false;
	    break;
	  }
	  char *tp = strchr(buf + 1, '\t');
	  if(!tp) {
	    fprintf(stderr,"Error in dbSNP file - bad chromosome header\n");
	    ok = false;
	    break;
	  }
	  *tp = 0;
	  char *tp1;
	  cbin = (int)strtoul(tp + 1, &tp1, 10);
	  if(*tp1 != '\t') {
	    fprintf(stderr,"Error in dbSNP file - bad chromosome header\n");
	    ok = false;
	    break;
	  }
	  max_bin = (int)strtoul(tp1 + 1, &tp, 10);
	  ctg_name = strdup(buf + 1);
	  HASH_FIND(hh, args.dbSNP, ctg_name, strlen(ctg_name), ctg);
	  if(ctg != NULL) {
	    fprintf(stderr,"Error in dbSNP file - duplicate contigs (%s)\n", ctg_name);
	    ok = false;
	    break;
	  }
	  ctg = malloc(sizeof(dbsnp_ctg));
	  ctg->name = ctg_name;
	  ctg->min_bin = cbin;
	  ctg->max_bin = max_bin;
	  ctg->bins = malloc(sizeof(dbsnp_bin) * (max_bin - cbin + 1));
	  HASH_ADD_KEYPTR(hh, args.dbSNP, ctg->name, strlen(ctg_name), ctg);
	  n_entries = name_buf_ptr = 0;
	  ix = -1;
	  n_ctgs++;
	  bins = ctg->bins;
	} else {
	  if(ctg == NULL) {
	    fprintf(stderr,"Error in dbSNP file - missing contig header\n");
	    ok = false;
	    break;
	  }
	  if(buf[0] == '+') {
	    if(n_entries) {
	      store_dbsnp_entries(bins, n_entries, name_buf_ptr, entries, name_buf);
	      n_bins++;
	      n_snps += n_entries;
	    }
	    char *tp;
	    int d = (int)strtoul(buf + 1, &tp, 10);
	    if(!d) d = 1;
	    cbin += d;
	    bins += d;
	    if(cbin > max_bin) {
	      fprintf(stderr,"Error in dbSNP file - too many bins for chromosome\n");
	      ok = false;
	      break;
	    }
	    n_entries = name_buf_ptr = 0;
	    ix = -1;
	  } else {
	    if(n_entries == 64) {
	      fprintf(stderr,"Error in dbSNP file - too many entries for bin (max 64)\n");
	      ok = false;
	      break;
	    }
	    if(l < 3 || l > 256) {
	      fprintf(stderr,"Error in dbSNP file - bad line length: %s\n", buf + 1);
	      ok = false;
	      break;
	    }
	    char tmp = buf[2];
	    buf[2] = 0;
	    char *tp;
	    int ix1 = (int)strtoul(buf, &tp, 16);
	    buf[2] = tmp;
	    int prefix_ix = ix1 >> 6;
	    if(prefix_ix > args.n_dbSNP_prefixes) {
	      fprintf(stderr,"Error in dbSNP file - invalid prefix\n");
	      ok = false;
	      break;
	    }
	    ix1 &= 63;
	    if(ix1 <= ix) {
	      fprintf(stderr,"Error in dbSNP file - entries out of order or invalid\n");
	      ok = false;
	      break;
	    }
	    ix = ix1;
	    int kx = 2;
	    if(prefix_ix == 0) {
	      if(l < 7) {
		ok = false;
		fprintf(stderr,"Error in dbSNP file - bad line length: %s\n", buf);
		break;
	      }
	      tmp = buf[4];
	      buf[4] = 0;
	      uint8_t ix_high = (int)strtoul(buf + 2, &tp, 16);
	      buf[4] = tmp;
	      tmp = buf[6];
	      buf[6] = 0;
	      uint8_t ix_low = (int)strtoul(buf + 4, &tp, 16);
	      buf[6] = tmp;
	      name_buf[name_buf_ptr++] = ix_high;
	      name_buf[name_buf_ptr++] = ix_low;
	      kx = 6;
	    }
	    int k = l - kx;
	    entries[n_entries++] = (k << 8) | (prefix_ix << 6) | ix;
	    uint8_t *tip = (uint8_t *)(buf + kx);
	    for(int j = 0; j < k; j++) name_buf[name_buf_ptr++] =  db_tab[(int)tip[j]];	    
	  }
	}
      }
      l = getline(&buf, &buf_size, fp);      
      if(l > 0 && buf[l - 1] == '\n') buf[--l] = 0;
    }
    if(n_entries) {
      store_dbsnp_entries(bins, n_entries, name_buf_ptr, entries, name_buf);
      n_bins++;
      n_snps += n_entries;
    }
  } else ok = false;
  fclose(fp);
  if(buf) free(buf);
  if(filter) {
    int i;
    while(waitpid(-1, &i, WNOHANG) > 0);
  }
  if(ok) fprintf(stderr,"Completed loading dbSNP (no. contigs %d, no. bins %d, no. SNPs %d\n", n_ctgs, n_bins, n_snps);
  else fprintf(stderr,"Error loading dbSNP\n");  
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out __unused__)
{
  args.hdr  = in;

  static struct option loptions[] = {
    {"output", required_argument, 0, 'o'},
    {"snps", required_argument, 0, 's'},
    {"dbsnp", required_argument, 0, 'D'},
    {"gzip", no_argument, 0, 'z'},
    {"bzip2", no_argument, 0, 'j'},
    {"xz", no_argument, 0, 'x'},
    {0,0,0,0}
  };
  int c;
  bool mult_comp = false;
  while ((c = getopt_long(argc, argv, "?ho:s:D:zjx",loptions,NULL)) >= 0) {
    switch (c) {
    case 'o':
      args.output_filename = optarg;
      break;
    case 's':
      args.snp_filename = optarg;
      break;
    case 'D':
      args.dbSNP_name = optarg;
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
  if (optind != argc) error(usage());
  init_files(&args);
  int ns = bcf_hdr_nsamples(args.hdr);
  assert(ns > 0);
  args.gt = malloc(sizeof(int) * ns);
  for(int i = 0; i < in->n[BCF_DT_ID]; i++) {
    if(!strcmp("PASS" , in->id[BCF_DT_ID][i].key)) args.pass_index = i;
    else if(!strcmp("GT" , in->id[BCF_DT_ID][i].key)) args.gt_index = i;
  }
  if(args.snp_filename != NULL) read_snp_file();
  if(args.dbSNP_name != NULL) read_dbSNP_file();
  FILE *fp = args.output_file;
  fputs("Chrom\tPos\tId", fp);
  for(int i = 0; i < ns; i++) {
    fprintf(fp, "\t%s", in->id[BCF_DT_SAMPLE][i].key);
  }
  fputc('\n', fp);
  return 1;
}

static const int base_tab[256] = {
  ['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4
};

static dbsnp_ctg *dbSNP_ctg;
static int curr_rid = -1;
static char dtab[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 0, 0, 0, 0, 0, 0 };

bcf1_t *process(bcf1_t *rec)
{
  int ns = bcf_hdr_nsamples(args.hdr);
  bcf_unpack(rec, BCF_UN_ALL);
  char *id = rec->d.id;
  char rs[512];
  void *dat_p = NULL;
  int dat_n = 0;
  FILE *fp = args.output_file;
  if(id[1] == 0 && id[0] == '.') {
    if(curr_rid != rec->rid) {
      const char *ctg = args.hdr->id[BCF_DT_CTG][rec->rid].key;
      if(args.dbSNP != NULL) {
	HASH_FIND(hh, args.dbSNP, ctg, strlen(ctg), dbSNP_ctg);
      }
      curr_rid = rec->rid;
    }
    if(dbSNP_ctg != NULL) {
      int x = rec->pos + 1;
      int bn = x >> 6;
      if(bn >= dbSNP_ctg->min_bin && bn <= dbSNP_ctg->max_bin) {
	dbsnp_bin *b = dbSNP_ctg->bins + bn - dbSNP_ctg->min_bin;
	int ix = x & 63;
	uint64_t mk = (uint64_t)1 << ix;
	if(b->mask & mk) {
	  uint64_t mk1 = b->mask & (mk - (uint64_t)1);
	  int i = 0, j = 0;
	  while(mk1) {
	    if(mk1 & (uint64_t)1) {
	      uint16_t en = b->entries[i++];
	      j += en >> 8;
	      if(!((en >> 6) & 3)) j += 2;
	    }
	    mk1 >>= 1;
	  }
	  char *tp = rs;
	  int prefix_id = (b->entries[i] >> 6) & 3;
	  unsigned char *tp1 = b->name_buf + j;
	  if((prefix_id--) == 0) {
	    prefix_id = (tp1[0] << 8) | tp1[1];
	    tp1+=2;
	  }
	  char *tp2 = args.dbSNP_prefix[prefix_id];
	  while(*tp2) *tp++ = *tp2++;
	  j = b->entries[i] >> 8;
	  for(int k = 0; k < j; k++) {
	    unsigned char z = *tp1++;
	    *tp++ = dtab[z >> 4];
	    *tp++ = dtab[z & 15];
	  }
	  *tp = 0;
	  id = rs;
	  //	  fprintf(stderr,"%s\n", rs);
	}
      }
    }    
  }
  if(id != NULL && (id[0] != '.' || id[1] != 0)) {
    bool passed = true;
    if(args.snp_hash) {
      snp *s;
      HASH_FIND_STR(args.snp_hash, id, s);
      if(!s) passed = false;
    }
    if(passed) {
      //      passed = false;
      passed = true;
      // Check overall filter
      for(int i = 0; i < rec->d.n_flt; i++) {
	if(rec->d.flt[i] == args.pass_index) {
	  passed = true;
	  break;
	}
      }
      int n_all = rec->n_allele;
      if(passed) {
	// Check alleles (only allow SNPs)
	if(n_all > 4) passed = false;
	else {
	  for(int i = 0; i < rec->n_allele; i++) {
	    char *p = rec->d.allele[i];
	    if(p[1] || !base_tab[(int)p[0]]) {
	      passed = false;
	      break;
	    }
	  }
	}
      }
      if(passed) {
	// Get filter tag
	int ne = bcf_get_format_values(args.hdr, rec, "FT", &dat_p, &dat_n, BCF_HT_STR);
	// Find GT tag
	int gt_i = -1;
	bcf_fmt_t *fmt = rec->d.fmt;
	for(int i = 0; i < (int)rec->n_fmt; i++) {
	  if(!fmt[i].p) continue;
	  if(fmt[i].id == args.gt_index) {
	    gt_i = i;
	    break;
	  }
	}
	if(gt_i >= 0) {
	  int sz = ne / ns;
	  char *flt = dat_p;
	  bcf_fmt_t *fmt = rec->d.fmt + gt_i;
	  passed = false;
	  for(int i = 0; i < ns; i++) {
	    args.gt[i] = 0;
	    switch(fmt->type) {
	      case BCF_BT_INT8:
		{
		  if(fmt->n == 2) {
		    int8_t *p = (int8_t *)(fmt->p + i * fmt->size);
		    if(p[0] != bcf_int8_vector_end && p[1] != bcf_int8_vector_end) {
		      int a1 = p[0] >> 1;
		      int a2 = p[1] >> 1;
		      if(a1 < 1 || a1 > n_all || a2 < 1 || a2 > n_all) args.gt[i] = 0;
		      else args.gt[i] = (a1 << 4) | a2;
		    }
		  }
		}
		break;
	      case BCF_BT_INT16:
		{
		  if(fmt->n == 2) {
		    int16_t *p = (int16_t *)(fmt->p + i * fmt->size);
		    if(p[0] != bcf_int16_vector_end && p[1] != bcf_int16_vector_end) {
		      int a1 = p[0] >> 1;
		      int a2 = p[1] >> 1;
		      if(a1 < 1 || a1 > n_all || a2 < 1 || a2 > n_all) args.gt[i] = 0;
		      else args.gt[i] = (a1 << 4) | a2;
		    }
		  }
		}
		break;
	      case BCF_BT_INT32:
		{
		  if(fmt->n == 2) {
		    int32_t *p = (int32_t *)(fmt->p + i * fmt->size);
		    if(p[0] != bcf_int32_vector_end && p[1] != bcf_int32_vector_end) {
		      int a1 = p[0] >> 1;
		      int a2 = p[1] >> 1;
		      if(a1 < 1 || a1 > n_all || a2 < 1 || a2 > n_all) args.gt[i] = 0;
		      else args.gt[i] = (a1 << 4) | a2;
		    }
		  }
		}
		break;
	    }
	    if(flt != NULL && strcmp("PASS", flt + i * sz))  args.gt[i] = 0;
	    if(args.gt[i]) passed = true;
	  }
	  if(passed) {
	    fprintf(fp, "%s\t%d\t%s", args.hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + 1, id);
	    for(int i = 0; i < ns; i++) {
	      const int gt = args.gt[i];
	      if(gt > 0) fprintf(fp, "\t%s%s", rec->d.allele[(gt >> 4) - 1], rec->d.allele[(gt & 7) - 1]);
	      else fputs("\t00", fp);
	    }
	    fputc('\n',fp);
	  }
	}
      }
    }
  }
  return NULL;
}

void destroy(void)
{
  close_files(&args);
}
