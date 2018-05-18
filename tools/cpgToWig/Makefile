#==================================================================================================
# PROJECT: gemBS
# FILE: Makefile
# DATE: 01/11/2015
# AUTHOR(S): Marcos Fernandez <mfernandez@cnag.crg.eu>
# DESCRIPTION: Builds the external tools for gemBS
#==================================================================================================


# Definitions

CC=gcc

ROOT_PATH=..

TOOLS=cpgToWig
FOLDER_BIN=../bin/
TOOLS_SRC=$(addsuffix .c, $(TOOLS))
TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(TOOLS))
GEN_LIBS:=-I../gen/include -L../gen/libsrc
LIBS:=-lgen -lm


ifeq ($(HAVE_ZLIB),1)
LIBS:=$(LIBS) -lz
endif
ifeq ($(HAVE_BZLIB),1)
LIBS:=$(LIBS) -lbz2
endif

all: TOOLS_FLAGS=-O2 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: $(TOOLS_BIN)
	
static: TOOLS_FLAGS=-O2 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS) -static
static: $(TOOLS_BIN)

debug: TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(TOOLS_BIN)

$(TOOLS_BIN): $(TOOLS_SRC)
	$(CC) --std=gnu11  $(TOOLS_FLAGS) -o $@ $(notdir $@).c $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(GEN_LIBS) $(LIBS) $(EXTRA_LIBS)
