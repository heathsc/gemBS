#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Builds all the sources of the GT-library, and the library itself (.a)
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

MODULES= gem_tools gt_commons gt_error gt_mm gt_ihash gt_shash gt_vector gt_string \
gt_input_file gt_pipe_io gt_fm

SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))
GT_LIB=$(FOLDER_LIB)/libgemtools.a

all: GEM_TOOLS_FLAGS=$(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZATION_FLAGS) $(ARCH_FLAGS_OPTIMIZATION_FLAGS)
all: $(GT_LIB)

debug: GEM_TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(GT_LIB)

$(GT_LIB): $(OBJS)
	$(AR) -rcs $(GT_LIB) $(FOLDER_BUILD)/*.o

$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_TOOLS_FLAGS) $(INCLUDE_FLAGS) $< -c -o $@
