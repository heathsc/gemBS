#ifndef BS_CALL_H

#include <assert.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

#include "dbSNP.h"
#define BS_CALL_VERSION "2.1.7"

#define STRING_EXP(tok) #tok
#define STRING(tok) STRING_EXP(tok)

#define DEFAULT_MAPQ_THRESH 20
#define DEFAULT_MAX_TEMPLATE_LEN 1000
#define DEFAULT_UNDER_CONVERSION 0.01
#define DEFAULT_OVER_CONVERSION 0.05
#define DEFAULT_REF_BIAS 2

#define MAX_GAUSS_N 256

// IMPORTANT!
//
// MIN_QUAL must not be lower than 1, and MAX_QUAL must be less than FLT_QUAL
// and FLT_QUAL must be less than 64.
//
#define MAX_QUAL 43
#define MIN_QUAL 20

#define FLT_QUAL 63

// #define QUAL_CONV 33
#define MAX_ITER 15
#define ITER_FIN (1.0e-8)

#define LOG10 (2.30258509299404568402)
#define LOG2 (0.69314718055994530942)

#define LFACT_STORE_SIZE 256

#define GET_QUAL(x) ((x) >> 2)
#define GET_BASE(x) ((x) & 3)

typedef enum {stats_all = 0, stats_passed} stats_cat;
typedef enum {all_sites = 0, variant_sites, CpG_ref_sites, CpG_nonref_sites} qual_cat;
typedef enum {mut_AC = 0, mut_AG, mut_AT, mut_CA, mut_CG, mut_CT, mut_GA, mut_GC, mut_GT, mut_TA, mut_TC, mut_TG, mut_no} stats_mut;
typedef enum {graphical_model, maximum_likelihood} BS_Caller;
typedef enum {base_none = 0, base_trim, base_clip, base_overlap, base_lowqual} base_filter_types;
typedef enum { NON_CONVERTED, STRAND_C2T, STRAND_G2A } gt_bs_strand;
typedef enum { gt_flt_none = 0, gt_flt_unmapped, gt_flt_qc, gt_flt_secondary, gt_flt_mate_unmapped, gt_flt_duplicate, gt_flt_nopos, gt_flt_nomatepos, gt_flt_mismatch_chr, gt_flt_orientation, gt_flt_insert_size, gt_flt_noseq, gt_flt_mapq,gt_flt_not_correctly_aligned } gt_filter_reason;
typedef enum { FORWARD, REVERSE, UNKNOWN } gt_strand;
typedef enum { NONE, GZIP, BZIP2, BGZIP } gt_compression;

typedef enum { MISMS, INS, DEL, SOFT } gt_misms_t;
typedef struct {
  gt_misms_t misms_type;
  uint32_t position;
  union {
    uint32_t size;
    char base;
  };
} gt_misms;

typedef struct {
  uint32_t forward_position;
  uint32_t reverse_position;
  uint32_t reference_span[2];
  gt_vector* read[2];
  gt_vector* mismatches[2];
  uint8_t mapq[2];
  gt_strand orientation;
  gt_bs_strand bs_strand;
} align_details;

typedef struct {
	uint64_t snps[2];
	uint64_t indels[2];
	uint64_t multi[2];
	uint64_t dbSNP_sites[2]; // How many dbSNP sites are covered
	uint64_t dbSNP_var[2]; // How many dbSNP sites are variant
	uint64_t CpG_ref[2];
	uint64_t CpG_nonref[2];
	int nbins;
	uint8_t *gc;
} gt_ctg_stats;

typedef struct {
	uint64_t coverage;
	uint64_t var;
	uint64_t CpG[2];
	uint64_t CpG_inf[2];
	uint64_t all;
	uint64_t gc_pcent[101];
	UT_hash_handle hh;
} gt_cov_stats;

typedef struct {
	char *name;
	struct _region *curr_reg;
	int bam_tid;
	int fai_id;
	int vcf_rid;
	uint32_t seq_len;
	uint32_t start_pos;
	uint32_t end_pos;
	uint16_t *seq;
	gt_ctg_stats *ctg_stats;
} ctg_t;

typedef struct _region {
	ctg_t *ctg;
	uint32_t start;
	uint32_t stop;
} region_t;

typedef struct {
	uint64_t conv_cts[4];
} meth_cts;

typedef struct {
	uint64_t cts[2];
} fstats_cts;

typedef struct {
	uint64_t snps[2];
	uint64_t indels[2];
	uint64_t multi[2];
	uint64_t dbSNP_sites[2]; // How many dbSNP sites are covered
	uint64_t dbSNP_var[2]; // How many dbSNP sites are variant
	uint64_t CpG_ref[2];
	uint64_t CpG_nonref[2];
	uint64_t mut_counts[12][2];
	uint64_t dbSNP_mut_counts[12][2];
	uint64_t qual[4][256];
	uint64_t filter_cts[15];
	uint64_t filter_bases[15];
	uint64_t base_filter[5];
	gt_vector *meth_profile;
	gt_vector *fs_stats;
	gt_vector *qd_stats;
	gt_vector *mq_stats;
	uint64_t filter_counts[2][32];
	double CpG_ref_meth[2][101];
	double CpG_nonref_meth[2][101];
	gt_cov_stats *cov_stats;
} bs_stats;

typedef struct {
	double e, k, ln_k, ln_k_half, ln_k_one;
} qual_prob;

typedef struct {
	uint64_t counts[8];
	int qual[8]; // Average quality per base type
	double gt_prob[10]; // Genotype log probabilities (Log10)
	double fisher_strand; // Allele strand bias LR (Log10)
	int mq; // RMS Mapping quality
	int aq; // Average base quality
	uint8_t max_gt;
} gt_meth;

typedef struct {
	gt_meth gtm;
	bool ready;
	bool skip;
} gt_vcf;

typedef struct _base_counts {
	struct _base_counts *next;
	int8_t idx[MAX_QUAL + 1];
	uint32_t counts[MAX_QUAL + 1];
} base_counts;

typedef struct {
	uint32_t counts[2][8];
	uint32_t n;
	float quality[8];
	float mapq2;
//	uint8_t *seq;
//	size_t seq_size;
//	size_t seq_idx;
} pileup;

typedef struct {
	align_details *al;
	uint32_t alignment_flag;
	uint32_t ix;
	gt_string *tag;
	UT_hash_handle hh;
} align_hash;

#define VCF_FLT_PASS 0
#define VCF_FLT_FAIL 1
#define VCF_FLT_MAC1 2
#define VCF_INFO_CX 3
#define VCF_FMT_GT 4
#define VCF_FMT_FT 5
#define VCF_FMT_GL 6
#define VCF_FMT_GQ 7
#define VCF_FMT_DP 8
#define VCF_FMT_MQ 9
#define VCF_FMT_QD 10
#define VCF_FMT_MC8 11
#define VCF_FMT_AMQ 12
#define VCF_FMT_CS 13
#define VCF_FMT_CG 14
#define VCF_FMT_FS 15
#define VCF_FMT_CX VCF_INFO_CX

typedef struct {
	uint32_t x;
	uint32_t max_pos;
	align_details *al;
	gt_vector *orig_pos[2];
} mprof_thread_t;

typedef struct {
	pthread_t thr;
	struct _sr_param *param;
	gt_meth *gtm;
	pileup *cts;
	uint32_t x;
	uint32_t y;
	int step;
	bool ready;
} cthread_par;

#define N_MPROF_BUFFERS 256

typedef struct {
	uint32_t y_waiting;
	ctg_t *ctg_waiting;
	gt_vector *align_list_waiting;
	gt_vector *free_list_waiting;
	bool process_end;
	bool print_end;
	bool mprof_end;
	bool calc_end;
	gt_vcf *vcf;
	uint32_t vcf_x;
	ctg_t *vcf_ctg;
	int vcf_size, vcf_n;
	pthread_mutex_t print_mutex;
	pthread_cond_t print_cond1;
	pthread_cond_t print_cond2;
	pthread_mutex_t vcf_mutex;
	pthread_cond_t vcf_cond;
	pthread_mutex_t process_mutex;
	pthread_cond_t process_cond1;
	pthread_cond_t process_cond2;
	pthread_mutex_t mprof_mutex;
	pthread_cond_t mprof_cond1;
	pthread_cond_t mprof_cond2;
	pthread_mutex_t calc_mutex;
	pthread_cond_t calc_cond1;
	pthread_cond_t calc_cond2;
	uint32_t n_contigs;
	uint32_t n_regions;
	faidx_t *seq_idx;
	ctg_t **contigs;
	int *tid2id;
	region_t *regions;
	region_t *curr_region;
	gt_string *ref;
	gt_string *ref1;
	dbsnp_header_t *dbSNP_hdr;
	bs_stats *stats;
	htsFile *vcf_file;
	bcf_hdr_t *vcf_hdr;
	int vcf_ids[17];
	FILE *json_file;
	bam_hdr_t *sam_header;
	htsFile *sam_file;
	hts_idx_t *sam_idx;
	cthread_par *calc_threads;
	int n_calc_threads;
	int calc_threads_complete;
	uint8_t flt_tab[768];
	mprof_thread_t mprof_thread[N_MPROF_BUFFERS];
	int mprof_read_idx;
	int mprof_write_idx;
} work_t;

typedef struct {
	char *flt_name[4];
	bool gt_het[10];
	double logp[100];
	double lfact_store[LFACT_STORE_SIZE];
} defs_t;

#define CALC_THREADS 0
#define INPUT_THREADS 1
#define OUTPUT_THREADS 2

typedef struct _sr_param {
	char *input_file;
	char *name_reference_file;
	char *output_file;
	char *sample_name;
	char *dbSNP_name;
	char *report_file;
	char *contig_bed;
	char *contig_sizes;
	bool mmap_input;
	/* Control flags */
	bool keep_duplicates;
	bool ignore_duplicates;
	bool keep_unmatched;
	bool haploid;
	bool verbose;
	bool blank_trim;
	bool all_positions;
	bool explicit_thread_distribution;
	bool benchmark_mode;
	uint8_t out_file_type;
	BS_Caller caller;
	uint32_t left_trim[2];
	uint32_t right_trim[2];
	int num_threads[3];
	uint8_t mapq_thresh;
	uint8_t min_qual;
	uint32_t max_template_len;
	double under_conv, over_conv;
	double ref_bias;

	// Definitions and tables
	defs_t defs;

	// Shared workspace
	work_t work;

} sr_param;


#define lfact2(x,lfs) ((x) < (LFACT_STORE_SIZE) ? lfs[x] : lgamma((double)((x) + 1)))

void fill_base_prob_table(void);
void calc_gt_prob(gt_meth *gt, const sr_param * const param, char rf);
void lfact_store_init(double * const lfact_store);
double fisher(int * const c, const double * const lfact_store);
void init_stats(sr_param *par);
void output_stats(sr_param *par);
void print_vcf_header(sr_param * const param, bam_hdr_t * hdr);
void flush_vcf_entries(bcf1_t * bcf, const sr_param * const par);
void print_vcf_entry(bcf1_t * bcf, ctg_t * const ctg, gt_meth *gtm, const char *rf,
		const uint32_t x, const uint32_t xstart, bool skip, sr_param * const par);
void init_param(sr_param * const par);
void trim_read(gt_vector * const sqread, int left_trim, int right_trim);
void trim_soft_clips(align_details * const al, bs_stats * const stats, uint32_t * const trim_left, uint32_t * const trim_right);
void handle_overlap(align_details * const al, bs_stats * const stats, uint32_t * const trim_left, uint32_t * const trim_right);
align_details *get_new_align_details(gt_vector * const free_list);
void clear_align_details(align_details * const al);
void make_align_hash(align_hash * const ah, gt_string * const tag, const uint32_t alignment_flag, const uint32_t ix, align_details * const al);
align_hash *make_new_align_hash(gt_vector * const free_hash_list, gt_string * const tag, const uint32_t alignment_flag, const uint32_t ix, align_details * const al);
void add_al_hash_to_free_list(align_hash * const ah, gt_vector * const free_hash_list);
bool check_for_tag_in_hash(const align_hash * const hash, gt_string * const tag);
void meth_profile(const align_details * const al, const uint32_t x, gt_vector * const orig_pos[2], const int max_pos, sr_param * const par);
void call_genotypes_ML(ctg_t * const ctg, gt_vector * const align_list, const uint32_t x, const uint32_t y, sr_param * const param);
void init_calc_threads(sr_param * const param);
void join_calc_threads(sr_param * const param);
bool get_sequence_index(sr_param * const par);
bool get_sequence_string(ctg_t * const ctg, const uint32_t x, const uint32_t sz, ctg_t * const prev_ctg, gt_string * const ref, sr_param * const par);
bool process_sam_header(sr_param * const par, bam_hdr_t * hdr);
gt_status open_input_file(sr_param * const param);
void close_input_file(sr_param * const param);

uint32_t get_al_qual(align_details *al);
gt_status parse_arguments(int argc, char **argv, sr_param *const par);
bool load_sequence(ctg_t * const contig, faidx_t *idx, const bool calc_gc);
void free_sequence(ctg_t * const contig);
gt_status bs_call_process(sr_param * const param);
gt_status read_input(htsFile *const sam_input, gt_vector * const align_list,sr_param *param);
gt_status process_template_vector(gt_vector *align_list, ctg_t * const ctg, uint32_t y, sr_param *param);
int get_next_align_details(htsFile * const sam_input, bam_hdr_t *hdr, hts_itr_t *itr, bam1_t *b, align_details * al, const uint64_t thresh,
		const uint64_t max_template_len, bool keep_unmatched, bool ignore_dup, bool *reverse, gt_filter_reason * const filtered,
		uint32_t * const align_length, uint32_t * const alignment_flag);

#define BS_CALL_H 1
#endif
