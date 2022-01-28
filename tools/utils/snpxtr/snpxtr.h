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

#include "dbSNP.h"

#define LOG10 2.30258509299404568402

void error(const char *format, ...) HTS_NORETURN;

typedef struct {
	bcf_hdr_t *hdr;
	bcf_srs_t *sr;
	char *snplistname;
	char *outfilename;
	char *dbSNPfilename;
	void *snp_hash;
	int *gt;
	dbsnp_header_t *dbSNP_hdr;
	htsFile *outfile;
	uint64_t *cumul_len;
	kstring_t out_string;
	int threads;
	int pass_idx;
	int gt_idx;
	bool compress;
	bool tabix;
	bool md5;
} sargs_t;

#define ks_output(fp, s) { \
	int r; \
	if((fp)->format.compression != no_compression) r = bgzf_write((fp)->fp.bgzf, (s)->s, (s)->l); \
	else r = hwrite((fp)->fp.hfile, (s)->s, (s)->l); \
	if(r != (s)->l) error("output error writing to %s\n", (fp)->fn ? (fp)->fn : "<NULL>"); \
} \

const char *usage(void);
void handle_command_line(int argc, char *argv[], sargs_t * const args);
void init_params(sargs_t *const args);
void read_snplist(sargs_t * const args);
void process_input(sargs_t * const args);
htsFile *open_ofile(char ** const name, bool compress, sargs_t * const a);
void calc_file_md5(char * const name);

#endif // MEXTR_H_
