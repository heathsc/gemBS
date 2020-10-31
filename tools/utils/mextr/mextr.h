#ifndef MEXTR_H_
#define MEXTR_H_

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

#include "bbi_structs.h"

#define LOG10 2.30258509299404568402

#define DEFAULT_UNDER_CONV 0.01
#define DEFAULT_OVER_CONV 0.05
#define DEFAULT_MAPQ_THRESH 20
#define DEFAULT_BQ_THRESH 20
#define DEFAULT_REF_BIAS 2
#define MAX_QUAL 43
#define DEFAULT_SELECT_THRESH 20

#define READ_BUF_SIZE 1024
#define REC_BUF_SIZE 2048

void error(const char *format, ...) HTS_NORETURN;

typedef struct {
	uint64_t n_sites;
	uint64_t n_sites_pass;
} stats_t;

typedef enum {FMT_FT, FMT_MC8, FMT_AMQ, FMT_CX, FMT_AQ, FMT_MQ, FMT_GQ, FMT_GOF, FMT_GL} fmt_tag;
typedef enum {CPGMODE_COMBINED, CPGMODE_SEPARATE} cpg_mode;
typedef enum {SELECT_HOM, SELECT_HET} select_mode;
typedef enum {BEDMETHYL_NONE = 0, BEDMETHYL_CPG, BEDMETHYL_CHG, BEDMETHYL_CHH} bedmethyl_type;

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

#define RJ_OUTPUT_CPG 1
#define RJ_OUTPUT_NONCPG 2
#define RJ_OUTPUT_BEDMETHYL 4
#define RJ_ALL (RJ_OUTPUT_CPG | RJ_OUTPUT_NONCPG | RJ_OUTPUT_BEDMETHYL)
#define REC_SKIP 64
#define REC_READY 128

#define N_FTAGS 6

typedef struct _rec {
	int32_t rid;
	hts_pos_t pos;
	char ref;
	char cx[5];
	uint8_t tasks;
	uint8_t max_common_gt;
	fmt_store_t tags[N_FTAGS];
	gt_meth *sample_gt;
} rec_t;

typedef struct {
	bcf1_t *buf[READ_BUF_SIZE];
	uint64_t idx[READ_BUF_SIZE];
	int read_pos;
	int write_pos;
	pthread_mutex_t mut;
	pthread_cond_t cond[2];
} bcf1_buffer_t;

typedef struct {
	rec_t *buf[REC_BUF_SIZE];
	uint64_t first_index;
	pthread_mutex_t mut;
} rec_buffer_t;

typedef struct {
	bbi_cblock_t *cblocks; // Buffer for holding output bbi data for compression
	uint32_t pos;
	int n_cblocks;
	pthread_mutex_t mut;
	pthread_cond_t cond[3];
	bool end_of_input;
} cblock_buffer_t;

typedef struct {
	bcf_hdr_t *hdr;
	char *cpgfilename;
	char *noncpgfilename;
	char *reportfilename;
	char *bedmethyl;
	char *bedmethylnames[3];
	char *bigbednames[3];
	char *bigwignames[2];
	char *bedmethyl_track_line;
	char *bedmethyl_desc;
	double *sample_Q[3];
	double *sample_Q1[3];
	cpg_prob *sample_cpg;
	htsFile *cpgfile;
	htsFile *noncpgfile;
	FILE *reportfile;
	htsFile *bedmethylfiles[3];
	FILE *bigbedfiles[3];
	FILE *bigwigfiles[2];
	bbi_header_t *bbi_hdr[2];
	bbi_ctg_data_t *ctg_data;
	uint64_t *cumul_len;
	int *id_trans; // translate from contig ids in BCF/VCF to contig ids in bbi files
	kstring_t pr_str[3];
	bbi_global_data_t bb_global[5];
	stats_t *stats;
	cpg_mode mode;
	select_mode sel_mode;
	bcf_srs_t *sr;
	bcf1_buffer_t read_buf;
	rec_buffer_t rec_buf;
	cblock_buffer_t cblock_buf;
	int sel_thresh;
	int threads;
	int compress_threads;
	bool compress;
	bool common_gt;
	bool output_noncpg;
	bool header;
	bool proc_finished;
	bool input_finished;
	bool rec_finished;
	bool strand_specific;
	bool calc_md5;
	bool tabix;
	uint8_t job_mask;
	double min_prop;
	int min_num;
	int min_inform;
	int min_nc;
	double ref_bias;
	double under_conv;
	double over_conv;
	int bq_thresh;
	int mq_thresh;
} args_t;

typedef struct {
	args_t * args;
	void (*output)(args_t *const, const rec_t * const);
	uint8_t job_flag;
} thr_info_t;

typedef struct {
	args_t * args;
	int thread_idx;
} gthr_info_t;

void calc_gt_prob(gt_meth *gt, args_t *args, char rf);
void calc_cpg_meth(args_t * const args, int ns, cpg_prob *cpg, gt_meth *g1, gt_meth *g2);
double get_meth(gt_meth *g, int idx);
void output_cpg(args_t *const args, rec_t ** const lrec, const int idx);
void output_noncpg(args_t *const args, const rec_t * const rec);
void output_bedmethyl(args_t *const args, const rec_t * const rec);
void fill_base_prob_table(void);
void print_headers(args_t *args);
int calc_phred(double z);
double *get_prob_dist(int ns, double *Q[]);
extern char trans_base[256];

#define ks_output(fp, s) { \
	int r; \
	if((fp)->format.compression != no_compression) r = bgzf_write((fp)->fp.bgzf, (s)->s, (s)->l); \
	else r = hwrite((fp)->fp.hfile, (s)->s, (s)->l); \
	if(r != (s)->l) error("output error writing to %s\n", (fp)->fn ? (fp)->fn : "<NULL>"); \
} \

const char *usage(void);
void handle_command_line(int argc, char *argv[], args_t * const args);
void init_params(args_t *const args);
void init_files(args_t *a);
void close_files(args_t *a);
void init_stats(args_t *a);
void write_stats(args_t *a);
void *unpack_bcf_thread(void *p);
void *read_thread(void * const p);
rec_t *rec_init(const int ns);
void *handle_rec_buf(void *p);
void *cpg_thread(void *p);
void *output_thread(void *p);
void calc_file_md5(char * const name);

#endif // MEXTR_H_
