#ifndef MEXTR_H_
#define MEXTR_H_

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <htslib/vcf.h>

#define COMP_GZIP (1 << COMPRESS_GZIP)
#define COMP_BZIP2 (1 << COMPRESS_BZIP2)
#define COMP_XZ (1 << COMPRESS_XZ)

#define LOG10 2.30258509299404568402

#define DEFAULT_UNDER_CONV 0.01
#define DEFAULT_OVER_CONV 0.05
#define DEFAULT_MAPQ_THRESH 20
#define DEFAULT_BQ_THRESH 20
#define DEFAULT_REF_BIAS 2
#define MAX_QUAL 43
#define DEFAULT_SELECT_THRESH 20

void error(const char *format, ...) HTS_NORETURN;

typedef struct {
  uint64_t n_sites;
  uint64_t n_sites_pass;
} stats_t;

typedef enum {FMT_FT, FMT_MC8, FMT_AMQ, FMT_CX, FMT_AQ, FMT_MQ, FMT_GQ, FMT_GOF, FMT_GL} fmt_tag;
typedef enum {CPGMODE_COMBINED, CPGMODE_SEPARATE} cpg_mode;
typedef enum {SELECT_HOM, SELECT_HET} select_mode;
typedef enum {BEDMETHYL_CPG, BEDMETHYL_CHG, BEDMETHYL_CHH, BEDMETHYL_NONE} bedmethyl_type;
	      
typedef struct {
  bcf_hdr_t *hdr;
  char *cpgfilename;
  char *noncpgfilename;
  char *reportfilename;
  char *wigfilename;
  char *bedmethyl;
  char *bedmethylnames[3];
  char *bedmethyl_track_line;
  char *bedmethyl_desc;
  FILE *cpgfile;
  FILE *noncpgfile;
  FILE *wigfile;
  FILE *reportfile;
  FILE *bedmethylfiles[3];
  stats_t *stats;
  cpg_mode mode;
  select_mode sel_mode;
  int sel_thresh;
  int compress;
  bool common_gt;
  bool output_noncpg;
  bool header;
  double min_prop;
  int min_num;
  int min_inform;
  int min_nc;
  double ref_bias;
  double under_conv;
  double over_conv;
  int bq_thresh;
  int mq_thresh;
  bool append_mode;
} args_t;

typedef struct {
  void *dat_p;
  int dat_n;
  int ne;
} fmt_store_t;

typedef struct {
  char *tag;
  int type;
  fmt_store_t st[2];
} fmt_field_t;

typedef struct {
  int32_t counts[8];
  int32_t aqual[8]; // Average base quality
  double gt_prob[10]; // Genotype log probabilities (Log10)
  double cmeth[3], gmeth[3];
  double sum;
  uint8_t max_gt;
  bool skip;
} gt_meth;

typedef struct {
  double prob_best;
  double prob_cg;
  uint8_t max_gt[2];
  double m;
} cpg_prob;
  
void calc_gt_prob(gt_meth *gt, args_t *args, char rf);
void calc_cpg_meth(args_t *args, int ns, cpg_prob *cpg, gt_meth *g1, gt_meth *g2);
double get_meth(gt_meth *g, int idx);
void output_cpg(args_t *args, bcf1_t *rec, fmt_field_t *tags, gt_meth *sample_gt[], int idx, cpg_prob *sample_cpg, double *Q[]);
void output_bedmethyl(args_t *args, bcf1_t *rec, fmt_field_t *tags, gt_meth *sample_gt[], int idx);
void fill_base_prob_table(void);

#endif // MEXTR_H_
