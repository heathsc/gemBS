#ifndef BS_CALL_OPTIONS_H

/*
* bs_call options
*/
gt_option bs_call_options[] = {
  /* Operations */
  { '1', "haploid", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 1, true, "", "Assume genome is haploid"},
  { 'd', "keep-duplicates", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 1, true, "", "Don't merge duplicate reads"},
  { 101, "ignore-duplicates", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 1, true, "", "Ignore duplicate flag from SAM/BAM files"},
  { 'k', "keep-unmatched", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 1, true, "", "Don't discard reads that don't form proper pairs"},
  { 'R', "right-trim", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "<int> | <int>,<int>", "Bases to trim from right of read pair for both reads or for read 1,read 2"},
  { 'L', "left-trim", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "<int> | <int>,<int>", "Bases to trim from left of read pair for both reads or for read 1,read 2"},
  { 'B', "blank-trim", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 1, true, "", "Don't use trimmed bases for genotype estimation"},
  { 'q', "mapq-threshold", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "<int> ","Set MAPQ threshold for selecting reads (default "STRING(DEFAULT_MAPQ_THRESH)")"},
  { 'Q', "bq-threshold", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "<int> ","Set base quality threshold for calling (default "STRING(MIN_QUAL)")"},
  { 'l', "max-template-length", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "<int> ","Set maximum template length for a pair (default "STRING(DEFAULT_MAX_TEMPLATE_LEN)")"},
//  { 'T', "realign-tolerance", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "<int> ","Tolerance for realignment positions (default "STRING(DEFAULT_REALIGN_TOL)")"},
  /* I/O */
  { 'O', "output-type", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<b|u|z|v|>" , "b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]" },
  { 202, "report-file", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<output JSON file name>" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<output file name>" , "" },
  { 'n', "sample", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<sample name>" , "SAMPLE" },
  { 'r', "reference", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<file> (MultiFASTA/FASTA)" , "" },
  { 'C', "contig-bed", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<file> (BED)" , "" },
  { 's', "contig-sizes", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<file>" , "" },
  { 'D', "dbsnp", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<file> (dbSNP processed file)" , "" },
  { 'A', "all-positions", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2, true, "" , "" },
  { 203, "benchmark-mode", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2, true, "" , "" },
	/* Model */
  { 'c', "conversion", GT_OPT_REQUIRED, GT_OPT_INT, 3, true, "<float>,<float>","Set under and over conversion rates (default "STRING(DEFAULT_UNDER_CONVERSION)","STRING(DEFAULT_OVER_CONVERSION)")"},	
  { 303, "reference-bias", GT_OPT_REQUIRED, GT_OPT_STRING, 3, true, "<float>","Set bias to reference homozygote (default "STRING(DEFAULT_REF_BIAS)")"},
	/* Misc */
  { 'v', "verbose", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4, true, "", ""},
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 4, true, "<int> | <int>,<int>,<int>", "Set additional threads for calculation, input and output"},
  { 'h', "help", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4, true, "", ""},
  {  0, 0, 0, 0, 0, false, "", ""}
};

char* bs_call_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Operations",
  /*  2 */ "I/O",
  /*  3 */ "Model",
  /*  4 */ "Misc",
};

#define BS_CALL_OPTIONS_H 1
#endif
