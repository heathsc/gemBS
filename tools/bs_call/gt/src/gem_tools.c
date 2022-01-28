/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gem_tools.h"

GT_INLINE uint64_t gt_options_get_num_options(const gt_option* const options) {
  uint64_t num_options = 0, i = 0;
  while (options[i++].option_id != 0) ++num_options;
  return num_options;
}
GT_INLINE struct option* gt_options_adaptor_getopt(const gt_option* const options) {
  const uint64_t num_options = gt_options_get_num_options(options);
  struct option* menu_options = gt_malloc(sizeof(struct option)*(1+num_options));
  // Adapt all the records
  uint64_t i = 0;
  for (i=0;i<=num_options;++i) {
    menu_options[i].name = options[i].long_option;
    menu_options[i].has_arg = options[i].option_type;
    menu_options[i].flag = 0;
    menu_options[i].val = options[i].option_id;
  }
  return menu_options;
}
GT_INLINE gt_string* gt_options_adaptor_getopt_short(const gt_option* const options) {
  const uint64_t num_options = gt_options_get_num_options(options);
  gt_string* const options_short = gt_string_new(2*num_options);
  // Adapt all the short options
  uint64_t i = 0;
  for (i=0;i<num_options;++i) {
    const char short_option = options[i].option_id;
    if (options[i].option_id<128 && gt_is_alphanumeric(short_option)) {
      gt_string_append_char(options_short,short_option);
      if (options[i].option_type==GT_OPT_REQUIRED || options[i].option_type==GT_OPT_OPTIONAL) {
        gt_string_append_char(options_short,COLON);
      }
    }
  }
  gt_string_append_eos(options_short);
  return options_short;
}
GT_INLINE void gt_options_fprint_menu(
    FILE* const stream,const gt_option* const options,char* groups[],
    const bool print_description,const bool print_inactive) {
  const uint64_t num_options = gt_options_get_num_options(options);
  int64_t i, last_group = -1;
  for (i=0;i<num_options;++i) {
    if (!print_inactive && !options[i].active) continue;
    // Print group (if not printed yet)
    if (last_group!=options[i].group_id) {
      fprintf(stream,"    [%s]\n",groups[options[i].group_id]);
      last_group=options[i].group_id;
    }
    // Print Long Option
    fprintf(stream,"      --%s",options[i].long_option);
    // Print Short Option (if it has)
    const char short_option = options[i].option_id;
    if (options[i].option_id<128 && gt_is_alphanumeric(short_option)) {
      fprintf(stream,"|-%c",short_option);
    }
    // Print extra command line syntax info
    fprintf(stream," %s",options[i].command_info);
    // Print description (@print_description)
    if (print_description && !gt_streq(options[i].description,"")) {
      fprintf(stream," %s",options[i].description);
    }
	  fputc('\n',stream);
  }
}
GT_INLINE void gt_options_fprint_json_menu(
    FILE* const stream,const gt_option* const options,char* groups[],
    const bool print_description,const bool print_inactive) {
  const uint64_t num_options = gt_options_get_num_options(options);
  int64_t i;
  bool at_least_one_printed = false;
  fprintf(stream,"{ \n"); // Begin JSON record
  fprintf(stream,"\"numOptions\": %"PRIu64",\n",num_options);
  fprintf(stream,"\"options\": [ \n");
  for (i=0;i<num_options;++i) {
    if (!print_inactive && !options[i].active) continue;
    if(at_least_one_printed){
      fprintf(stream,",\n");
    }
    at_least_one_printed = true;
    fprintf(stream,"\t{ \n");
    // Print ID/Short Option
    fprintf(stream,"\t  \"ID\": %d,\n",options[i].option_id);
    // Print Long Option
    fprintf(stream,"\t  \"longOption\": \"%s\",\n",options[i].long_option);
    // Print Short Option
    const char short_option = options[i].option_id;
    if (options[i].option_id<128 && gt_is_alphanumeric(short_option)) {
      fprintf(stream,"\t  \"shortOption\": \"%c\",\n",short_option);
    } else {
      fprintf(stream,"\t  \"shortOption\": null,\n");
    }
    // Group
    fprintf(stream,"\t  \"group\": \"%s\",\n",groups[options[i].group_id]);
    // Option Type
    switch (options[i].option_type) {
      case GT_OPT_NO_ARGUMENT: fprintf(stream,"\t  \"optionType\": \"noArgument\",\n"); break;
      case GT_OPT_REQUIRED: fprintf(stream,"\t  \"optionType\": \"required\",\n"); break;
      case GT_OPT_OPTIONAL: fprintf(stream,"\t  \"optionType\": \"optional\",\n"); break;
    }
    // Argument Type
    switch (options[i].argument_type) {
      case GT_OPT_NONE: fprintf(stream,"\t  \"argumentType\": null,\n"); break;
      case GT_OPT_INT: fprintf(stream,"\t  \"argumentType\": \"int\",\n"); break;
      case GT_OPT_FLOAT: fprintf(stream,"\t  \"argumentType\": \"float\",\n"); break;
      case GT_OPT_CHAR: fprintf(stream,"\t  \"argumentType\": \"char\",\n"); break;
      case GT_OPT_STRING: fprintf(stream,"\t  \"argumentType\": \"string\",\n"); break;
      case GT_OPT_BOOL: fprintf(stream,"\t  \"argumentType\": \"bool\",\n"); break;
    }
    // Print extra command line syntax info
    fprintf(stream,"\t  \"commandInfo\": \"%s\"",options[i].command_info);
    // Print description (@print_description)
    if (print_description && !gt_streq(options[i].description,"")) {
      fprintf(stream,",\n\t  \"description\": \"%s\"",options[i].description);
    }
    fprintf(stream,"\n\t}");
  }
  fprintf(stream,"\n    ]\n");
  fprintf(stream,"}\n");
}


