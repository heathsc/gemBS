#==================================================================================================
# PROJECT: gemBS
# FILE: Makefile
# DATE: 10/05/2018
# AUTHOR(S): Simon Heath simon.heath@gmail.com
# DESCRIPTION: Builds the ReadNameClean utility for gemBS
#==================================================================================================


# Definitions

CC=gcc

ROOT_PATH=..

TOOLS=readNameClean
FOLDER_BIN=../bin
TOOLS_SRC=$(addsuffix .c, $(TOOLS))
TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(TOOLS))
LIBS:= -lm


all: TOOLS_FLAGS=-O3 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: $(TOOLS_BIN)
	
static: TOOLS_FLAGS=-O3 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS) -static
static: $(TOOLS_BIN)

debug: TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(TOOLS_BIN)

$(TOOLS_BIN): $(TOOLS_SRC)
	$(CC) --std=c99  $(TOOLS_FLAGS) -o $@ $(notdir $@).c $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS) $(EXTRA_LIBS)
