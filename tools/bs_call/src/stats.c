/*
 * stats.c
 *
 *  Created on: Nov 24, 2019
 *      Author: heath
 */

#include <stdio.h>

#include "gem_tools.h"
#include "bs_call.h"

static int sort_cov_stats(gt_cov_stats *a, gt_cov_stats *b) {
	if(a->coverage < b->coverage) return -1;
	else if(a->coverage > b->coverage) return 1;
	return 0;
}

void output_stats(sr_param *par) {
	static char *mut_type[] = {"A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"};
	static char *filter_names[] = {
			"Passed", "Unmapped", "QC_Flags", "SecondaryAlignment", "MateUnmapped", "Duplicate", "NoPosition", "NoMatePosition", "MismatchContig", "BadOrientation", "LargeInsertSize", "NoSequence", "LowMAPQ", "NotCorrectlyAligned", "PairNotFound"
	};
	static char *base_filters[] = {
			"Passed", "Trimmed", "Clipped", "Overlapping", "LowQuality"
	};
	FILE *fp = par->work.json_file;
	bs_stats *stats = par->work.stats;
	fprintf(fp, "{\n\t\"source\": \"bs_call_v2.1, under_conversion=%g, over_conversion=%g, mapq_thresh=%d, bq_thread=%d\",\n",
			par->under_conv, par->over_conv, par->mapq_thresh, par->min_qual);
	time_t cl = time(0);
	struct tm *tt = localtime(&cl);
	fprintf(fp, "\t\"date\": \"%02d/%02d/%04d\",\n", tt->tm_mday, tt->tm_mon + 1, tt->tm_year + 1900);
	fputs("\t\"filterStats\": {\n\t\t\"ReadLevel\": {\n", fp);
	fprintf(fp, "\t\t\t\"%s\": {\n\t\t\t\t\"Reads\": %" PRIu64 ",\n\t\t\t\t\"Bases\": %" PRIu64 "\n\t\t\t}", filter_names[0], stats->filter_cts[0], stats->filter_bases[0]);
	for(int i = 1; i < 15; i++) {
		if(stats->filter_cts[i] > 0) {
			fprintf(fp, ",\n\t\t\t\"%s\": {\n\t\t\t\t\"Reads\": %" PRIu64 ",\n\t\t\t\t\"Bases\": %" PRIu64 "\n\t\t\t}", filter_names[i], stats->filter_cts[i], stats->filter_bases[i]);
		}
	}
	fputs("\n\t\t},\n\t\t\"BaseLevel\": {\n", fp);
	fprintf(fp, "\t\t\t\"%s\": %" PRIu64, base_filters[0], stats->base_filter[0]);
	for(int i = 1; i < 5; i++) {
		if(stats->base_filter[i] > 0) {
			fprintf(fp, ",\n\t\t\t\"%s\": %" PRIu64, base_filters[i], stats->base_filter[i]);
		}
	}
	fputs("\n\t\t}\n\t},\n\t\"totalStats\": {\n", fp);
	fprintf(fp, "\t\t\"SNPS\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->snps[stats_all], stats->snps[stats_passed]);
	fprintf(fp, "\t\t\"Indels\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->indels[stats_all], stats->indels[stats_passed]);
	fprintf(fp, "\t\t\"Multiallelic\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->multi[stats_all], stats->multi[stats_passed]);
	if(par->work.dbSNP_hdr != NULL) {
		fprintf(fp, "\t\t\"dbSNPSites\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->dbSNP_sites[stats_all], stats->dbSNP_sites[stats_passed]);
		fprintf(fp, "\t\t\"dbSNPVariantSites\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->dbSNP_var[stats_all], stats->dbSNP_var[stats_passed]);
	}
	fprintf(fp, "\t\t\"RefCpG\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->CpG_ref[stats_all], stats->CpG_ref[stats_passed]);
	fprintf(fp, "\t\t\"NonRefCpG\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->CpG_nonref[stats_all], stats->CpG_nonref[stats_passed]);
	fputs("\t\t\"QCDistributions\": {\n", fp);
	fputs("\t\t\t\"FisherStrand\": ", fp);
	char term = '{';
	for(int i = 0; i < stats->fs_stats->used; i++) {
		fstats_cts *c = gt_vector_get_elm(stats->fs_stats, i, fstats_cts);
		if(c->cts[1] > 0) {
			fprintf(fp,"%c\n\t\t\t\t\"%d\": %" PRIu64, term, i, c->cts[1]);
			term = ',';
		}
	}
	if(term == '{') fputc(term, fp);
	fputs("\n\t\t\t},\n", fp);
	fputs("\t\t\t\"QualityByDepth\": ", fp);
	term = '{';
	for(int i = 0; i < stats->qd_stats->used; i++) {
		fstats_cts *c = gt_vector_get_elm(stats->qd_stats, i, fstats_cts);
		if(c->cts[0] + c->cts[1] > 0) {
			fprintf(fp,"%c\n\t\t\t\t\"%d\": {\"NonVariant\": %" PRIu64 ", \"Variant\": %" PRIu64 "}", term, i, c->cts[0], c->cts[1]);
			term = ',';
		}
	}
	if(term == '{') fputc(term, fp);
	fputs("\n\t\t\t},\n", fp);
	fputs("\t\t\t\"RMSMappingQuality\": ", fp);
	term = '{';
	for(int i = 0; i < stats->mq_stats->used; i++) {
		fstats_cts *c = gt_vector_get_elm(stats->mq_stats, i, fstats_cts);
		if(c->cts[0] + c->cts[1] > 0) {
			fprintf(fp,"%c\n\t\t\t\t\"%d\": {\"NonVariant\": %" PRIu64 ", \"Variant\": %" PRIu64 "}", term, i, c->cts[0], c->cts[1]);
			term = ',';
		}
	}
	if(term == '{') fputc(term, fp);
	fputs("\n\t\t\t}\n\t\t},\t\t\"VCFFilterStats\": {\n", fp);
	fprintf(fp,"\t\t\t\"PASS\": {\"NonVariant\": %" PRIu64 ", \"Variant\": %" PRIu64 "}", stats->filter_counts[0][0], stats->filter_counts[1][0]);
	for(int i = 1; i < 16; i++) {
		fputs(",\n\t\t\t", fp);
		int k = i;
		int f_ix = 0;
		char tmp = '\"';
		while(k) {
			if(k & 1) {
				fprintf(fp, "%c%s", tmp, par->defs.flt_name[f_ix]);
				tmp = ',';
			}
			k >>= 1;
			f_ix++;
		}
		fprintf(fp,"\": {\"NonVariant\": %" PRIu64 ", \"Variant\": %" PRIu64 "}", stats->filter_counts[0][i], stats->filter_counts[1][i]);
	}
	fputs("\n\t\t},\n", fp);
	HASH_SORT(stats->cov_stats, sort_cov_stats);
	fputs("\t\t\"coverage\": {\n",fp);
	fputs("\t\t\t\"All\": ", fp);
	int ix = 0;
	term = '{';
	for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
		if(gcov->all != 0) {
			if(!(ix++)) {
				fprintf(fp, "%c\n\t\t\t\t", term);
				term = ',';
			} else fputs(", ", fp);
			fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->all);
			ix %= 12;
		}
	}
	fputs("\n\t\t\t},\n\t\t\t\"Variant\": ", fp);
	term = '{';
	ix = 0;
	for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
		if(gcov->var != 0) {
			if(!(ix++)) {
				fprintf(fp, "%c\n\t\t\t\t", term);
				term = ',';
			} else fputs(", ", fp);
			fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->var);
			ix %= 12;
		}
	}
	fputs("\n\t\t\t},\n\t\t\t\"RefCpG\": ", fp);
	term = '{';
	ix = 0;
	for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
		if(gcov->CpG[0] != 0) {
			if(!(ix++)) {
				fprintf(fp, "%c\n\t\t\t\t", term);
				term = ',';
			} else fputs(", ", fp);
			fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->CpG[0]);
			ix %= 12;
		}
	}
	fputs("\n\t\t\t},\n\t\t\t\"RefCpGInf\": ", fp);
	term = '{';
	ix = 0;
	for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
		if(gcov->CpG_inf[0] != 0) {
			if(!(ix++)) {
				fprintf(fp, "%c\n\t\t\t\t", term);
				term = ',';
			} else fputs(", ", fp);
			fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->CpG_inf[0]);
			ix %= 12;
		}
	}
	fputs("\n\t\t\t},\n\t\t\t\"NonRefCpG\": ", fp);
	term = '{';
	ix = 0;
	for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
		if(gcov->CpG[1] != 0) {
			if(!(ix++)) {
				fprintf(fp, "%c\n\t\t\t\t", term);
				term = ',';
			} else fputs(", ", fp);
			fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->CpG[1]);
			ix %= 12;
		}
	}
	fputs("\n\t\t\t},\n\t\t\t\"NonRefCpGInf\": ", fp);
	term = '{';
	ix = 0;
	for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
		if(gcov->CpG_inf[1] != 0) {
			if(!(ix++)) {
				fprintf(fp, "%c\n\t\t\t\t", term);
				term = ',';
			} else fputs(", ", fp);
			fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->CpG_inf[1]);
			ix %= 12;
		}
	}
	fputs("\n\t\t\t},\n\t\t\t\"GC\": ", fp);
	term = '{';
	for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
		if(!gcov->all) continue;
		fprintf(fp, "%c\n\t\t\t\t\"%" PRIu64 "\": [\n\t\t\t\t\t", term, gcov->coverage);
		term = ',';
		for(int i = 0; i < 100; i++) {
			fprintf(fp, "%" PRIu64 ",", gcov->gc_pcent[i]);
			if((i & 15) == 15) fputs("\n\t\t\t\t\t", fp);
			else fputc(' ', fp);
		}
		fprintf(fp, "%" PRIu64 "\n\t\t\t\t]", gcov->gc_pcent[100]);
	}
	fputs("\n\t\t\t}\n\t\t},\n\t\t\"quality\": {\n", fp);
	fputs("\t\t\t\"All\": [\n\t\t\t\t", fp);
	for(int i = 0; i < 255; i++) {
		fprintf(fp, "%" PRIu64 ", ", stats->qual[all_sites][i]);
		if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
	}
	fprintf(fp, "%" PRIu64 "\n\t\t\t],\n", stats->qual[all_sites][255]);
	fputs("\t\t\t\"Variant\": [\n\t\t\t\t", fp);
	for(int i = 0; i < 255; i++) {
		fprintf(fp, "%" PRIu64 ",", stats->qual[variant_sites][i]);
		if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
		else fputc(' ', fp);
	}
	fprintf(fp, "%" PRIu64 "\n\t\t\t],\n", stats->qual[variant_sites][255]);
	fputs("\t\t\t\"RefCpG\": [\n\t\t\t\t", fp);
	for(int i = 0; i < 255; i++) {
		fprintf(fp, "%" PRIu64 ",", stats->qual[CpG_ref_sites][i]);
		if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
		else fputc(' ', fp);
	}
	fprintf(fp, "%" PRIu64 "\n\t\t\t],\n", stats->qual[CpG_ref_sites][255]);
	fputs("\t\t\t\"NonRefCpG\": [\n\t\t\t\t", fp);
	for(int i = 0; i < 255; i++) {
		fprintf(fp, "%" PRIu64 ",", stats->qual[CpG_nonref_sites][i]);
		if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
		else fputc(' ', fp);
	}
	fprintf(fp, "%" PRIu64 "\n\t\t\t]\n", stats->qual[CpG_nonref_sites][255]);
	fputs("\t\t},\n\t\t\"mutations\": {\n", fp);
	for(int mut = 0; mut < 11; mut++) {
		fprintf(fp, "\t\t\t\"%s\": { \"All\": %" PRIu64 ", \"Passed\": %" PRIu64 ", \"dbSNPAll\": %" PRIu64 ", \"dbSNPPassed\": %" PRIu64 " },\n",
				mut_type[mut], stats->mut_counts[mut][stats_all], stats->mut_counts[mut][stats_passed],
				stats->dbSNP_mut_counts[mut][stats_all], stats->dbSNP_mut_counts[mut][stats_passed]);
	}
	fprintf(fp, "\t\t\t\"%s\": { \"All\": %" PRIu64 ", \"Passed\": %" PRIu64 ", \"dbSNPAll\": %" PRIu64 ", \"dbSNPPassed\": %" PRIu64 " }\n",
			mut_type[11], stats->mut_counts[11][stats_all], stats->mut_counts[11][stats_passed],
			stats->dbSNP_mut_counts[11][stats_all], stats->dbSNP_mut_counts[11][stats_passed]);
	fputs("\t\t},\n\t\t\"methylation\": {\n", fp);
	fputs("\t\t\t\"AllRefCpg\": [\n\t\t\t\t", fp);
	for(int i = 0; i < 100; i++) {
		fprintf(fp, "%.8g, ", stats->CpG_ref_meth[stats_all][i]);
		if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
	}
	fprintf(fp, "%.8g\n\t\t\t],\n", stats->CpG_ref_meth[stats_all][100]);
	fputs("\t\t\t\"PassedRefCpg\": [\n\t\t\t\t", fp);
	for(int i = 0; i < 100; i++) {
		fprintf(fp, "%.8g, ", stats->CpG_ref_meth[stats_passed][i]);
		if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
	}
	fprintf(fp, "%.8g\n\t\t\t],\n", stats->CpG_ref_meth[stats_passed][100]);
	fputs("\t\t\t\"AllNonRefCpg\": [\n\t\t\t\t", fp);
	for(int i = 0; i < 100; i++) {
		fprintf(fp, "%.8g, ", stats->CpG_nonref_meth[stats_all][i]);
		if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
	}
	fprintf(fp, "%.8g\n\t\t\t],\n", stats->CpG_nonref_meth[stats_all][100]);
	fputs("\t\t\t\"PassedNonRefCpg\": [\n\t\t\t\t", fp);
	for(int i = 0; i < 100; i++) {
		fprintf(fp, "%.8g, ", stats->CpG_nonref_meth[stats_passed][i]);
		if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
	}
	fprintf(fp, "%.8g\n\t\t\t]", stats->CpG_nonref_meth[stats_passed][100]);
	int nr = gt_vector_get_used(stats->meth_profile);
	if(nr) {
//		for(uint32_t k = 0; k < nr; k++) {
//			meth_cts *mc = gt_vector_get_mem(stats->meth_profile, meth_cts);
//			fprintf(stderr,"%u: %lu %lu %lu %lu\n", k, mc[k].conv_cts[0], mc[k].conv_cts[1], mc[k].conv_cts[2], mc[k].conv_cts[3]);
//		}
		fputs(",\n\t\t\t\"NonCpGreadProfile\": ", fp);
		term ='[';
		for(int i = 1; i < nr; i++) {
			meth_cts *mc = gt_vector_get_elm(stats->meth_profile, i, meth_cts);
			fprintf(fp, "%c\n\t\t\t\t[ %" PRIu64 ", %" PRIu64 ", %" PRIu64 ", %" PRIu64 " ]", term, mc->conv_cts[0], mc->conv_cts[1], mc->conv_cts[2], mc->conv_cts[3]);
			term = ',';
		}
		fputs("\n\t\t\t]", fp);
	}
	fputs("\n\t\t}\n\t},\n\t\"contigStats\": ", fp);
	term = '{';
	for(int i = 0; i < par->work.n_contigs; i++) {
		ctg_t * const ctg = par->work.contigs[i];
		gt_ctg_stats * const gs = ctg->ctg_stats;
		if(gs == NULL || gs->snps[stats_all] == 0) continue;
		fprintf(fp, "%c\n\t\t\"%s\": {\n", term, ctg->name);
		term = ',';
		fprintf(fp, "\t\t\t\"SNPS\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->snps[stats_all], gs->snps[stats_passed]);
		fprintf(fp, "\t\t\t\"Indels\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->indels[stats_all], gs->indels[stats_passed]);
		fprintf(fp, "\t\t\t\"Multiallelic\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->multi[stats_all], gs->multi[stats_passed]);
		if(par->work.dbSNP_hdr!= NULL) {
			fprintf(fp, "\t\t\t\"dbSNPSites\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->dbSNP_sites[stats_all], gs->dbSNP_sites[stats_passed]);
			fprintf(fp, "\t\t\t\"dbSNPVariantSites\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->dbSNP_var[stats_all], gs->dbSNP_var[stats_passed]);
		}
		fprintf(fp, "\t\t\t\"RefCpG\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->CpG_ref[stats_all], gs->CpG_ref[stats_passed]);
		fprintf(fp, "\t\t\t\"NonRefCpG\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t}\n\t\t}", gs->CpG_nonref[stats_all], gs->CpG_nonref[stats_passed]);
	}
	fputs("\n\t}\n}\n", fp);
}

void init_stats(sr_param *par) {
	par->work.json_file = fopen(par->report_file, "w");
	if(par->work.json_file == NULL) fprintf(stderr,"Could not open report file '%s' for output: %s\n", par->report_file, strerror(errno));
	else {
		bs_stats * const stats = par->work.stats = calloc((size_t)1, sizeof(bs_stats));
		stats->meth_profile = gt_vector_new(256, sizeof(meth_cts));
		stats->qd_stats = gt_vector_new(256, sizeof(fstats_cts));
		stats->fs_stats = gt_vector_new(256, sizeof(fstats_cts));
		stats->mq_stats = gt_vector_new(256, sizeof(fstats_cts));
		memset(stats->meth_profile->memory, 0, sizeof(meth_cts) * stats->meth_profile->elements_allocated);
		memset(stats->qd_stats->memory, 0, sizeof(fstats_cts) * stats->qd_stats->elements_allocated);
		memset(stats->fs_stats->memory, 0, sizeof(fstats_cts) * stats->fs_stats->elements_allocated);
		memset(stats->mq_stats->memory, 0, sizeof(fstats_cts) * stats->mq_stats->elements_allocated);
	}
}
