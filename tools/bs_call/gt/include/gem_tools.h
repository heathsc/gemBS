/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GEM_TOOLS_H_
#define GEM_TOOLS_H_

// Essentials
#include "gt_essentials.h"

/*
 * Options (Tools Menu)
 */
typedef enum { GT_OPT_NO_ARGUMENT=no_argument, GT_OPT_REQUIRED=required_argument, GT_OPT_OPTIONAL=optional_argument } gt_option_t;
typedef enum { GT_OPT_NONE, GT_OPT_INT, GT_OPT_FLOAT, GT_OPT_CHAR, GT_OPT_STRING, GT_OPT_BOOL } gt_option_argument_t;
typedef struct {
  int option_id;       // Integer ID or short character option
  char* long_option;   // Long string option
  gt_option_t option_type;            // Option type
  gt_option_argument_t argument_type; // Type of the argument
  uint64_t group_id;   // Label of the group it belongs to (zero if none)
  bool active;         // Enable/Disable option
  char* command_info;  // Extra command line syntax info
  char* description;   // Brief description
} gt_option;

GT_INLINE uint64_t gt_options_get_num_options(const gt_option* const options);
GT_INLINE struct option* gt_options_adaptor_getopt(const gt_option* const options);
GT_INLINE gt_string* gt_options_adaptor_getopt_short(const gt_option* const options);
GT_INLINE void gt_options_fprint_menu(
    FILE* const stream,const gt_option* const options,char* gt_filter_groups[],
    const bool print_description,const bool print_inactive);
GT_INLINE void gt_options_fprint_json_menu(
    FILE* const stream,const gt_option* const options,char* groups[],
    const bool print_description,const bool print_inactive);

#endif /* GEM_TOOLS_H_ */
