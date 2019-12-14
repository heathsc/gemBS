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

TOOLS=gemBS_cat readNameClean md5_fasta
FOLDER_BIN=../bin
TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(TOOLS))
LIBS:= -lm


all: TOOLS_FLAGS=-O3 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: $(TOOLS_BIN)
	
static: TOOLS_FLAGS=-O3 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS) -static
static: $(TOOLS_BIN)

debug: TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(TOOLS_BIN)

clean: 
	rm -f *~ *.o *.a

utils.o: utils.c utils.h
	$(CC) $(TOOLS_FLAGS) -c utils.c

$(FOLDER_BIN)/gemBS_cat: gemBS_cat.c utils.o
	$(CC) $(TOOLS_FLAGS) -o $@ gemBS_cat.c utils.o $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS) $(EXTRA_LIBS)

$(FOLDER_BIN)/md5_fasta: md5_fasta.c utils.o
	$(CC) $(TOOLS_FLAGS) -o $@ md5_fasta.c utils.o $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS) $(EXTRA_LIBS)-lcrypto

$(FOLDER_BIN)/readNameClean: readNameClean.c utils.o
	$(CC) $(TOOLS_FLAGS) -o $@ readNameClean.c utils.o $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS) $(EXTRA_LIBS)
