#============================================================================
# PROJECT: gemBS
# FILE: Makefile
# DATE: 01/01/2015
# AUTHOR(S): Marcos Fernandez Callejo <mfernandez@cnag.crg.eu>
# DESCRIPTION: Top level makefile
#============================================================================

# Definitions
ROOT_PATH=$(CURDIR)

# samtools and bcftools definitions
SAMTOOLS_VERSION=1.9
BCFTOOLS_VERSION=1.9
UCSCTOOLS_VERSION=v378
SAMTOOLS_DIR=samtools
BCFTOOLS_DIR=bcftools
UCSCTOOLS_DIR=userApps
SAMTOOLS=$(SAMTOOLS_DIR)/samtools
BCFTOOLS=$(BCFTOOLS_DIR)/bcftools
SAMTOOLS_TAR=samtools-$(SAMTOOLS_VERSION).tar.bz2
BCFTOOLS_TAR=bcftools-$(BCFTOOLS_VERSION).tar.bz2
UCSCTOOLS_TAR=userApps.$(UCSCTOOLS_VERSION).src.tgz
SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/$(SAMTOOLS_VERSION)/$(SAMTOOLS_TAR)
BCFTOOLS_URL=https://github.com/samtools/bcftools/releases/download/$(BCFTOOLS_VERSION)/$(BCFTOOLS_TAR)
UCSCTOOLS_URL=http://hgdownload.cse.ucsc.edu/admin/exe/$(UCSCTOOLS_TAR)

MACHTYPE:=$(shell uname -m)
ifneq (,$(findstring -,$(MACHTYPE)))
    MACHTYPE:=$(shell uname -m)
endif

UCSCTOOLS_LIB=$(UCSCTOOLS_DIR)/kent/src/lib/$(MACHTYPE)/jkweb.a
UCSCTOOLS_HTSLIB=$(UCSCTOOLS_DIR)/kent/src/htslib/libhts.a

#UNAME := $(shell uname -s)
#ifeq ($(UNAME),Linux)
#  ARCHIVE=http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v369
#endif
#ifeq ($(UNAME),Darwin)
# ARCHIVE=http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64
#endif

		  
FOLDER_BIN=bin

all: setup gem3 _samtools _bcftools kent _bs_call _utils

all_static: setup gem3-static _samtools _bcftools kent _bs_call _utils

all_debug: setup gem3-debug _samtools _bcftools kent _bs_call _utils

_samtools: $(SAMTOOLS_DIR) $(SAMTOOLS_DIR)/config.mk 
	$(MAKE) --directory=$(SAMTOOLS_DIR) all all-htslib

_bcftools: $(BCFTOOLS_DIR) $(BCFTOOLS_DIR)/config.h
	$(MAKE) $(FOLDER_BIN)/bcftools

kent: $(FOLDER_BIN)/wigToBigWig $(FOLDER_BIN)/bedToBigBed

_bs_call: bs_call/src/Makefile.mk
	$(MAKE) --directory=bs_call

_utils:
	$(MAKE) --directory=utils

setup:
	@mkdir -p $(FOLDER_BIN) 

bs_call/src/Makefile.mk: bs_call/src/Makefile.mk.in bs_call/gt/Makefile.mk.in bs_call/configure _samtools
	cd bs_call; ./configure ${BS_CALL_CONFIG} --with-htslib=../${SAMTOOLS_DIR}/htslib-${SAMTOOLS_VERSION}

gem3: gem3-mapper/Makefile.mk
	$(MAKE) --directory=gem3-mapper
	
gem3-static: gem3-mapper/Makefile.mk
	$(MAKE) --directory=gem3-mapper static
	
gem3-debug: gem3-mapper/Makefile.mk
	$(MAKE) --directory=gem3-mapper debug
	
gem3-mapper/Makefile.mk: gem3-mapper/Makefile.mk.in gem3-mapper/configure
	cd gem3-mapper; ./configure ${GEM3_CONFIG}

$(SAMTOOLS_DIR)/config.mk:
	cd $(SAMTOOLS_DIR); ./configure --disable-lzma

$(FOLDER_BIN)/bcftools: $(BCFTOOLS_DIR) $(BCFTOOLS_DIR)/plugins/mextr.c $(BCFTOOLS_DIR)/plugins/snpxtr.c
	$(MAKE) --directory=$(BCFTOOLS_DIR) all

$(BCFTOOLS_DIR)/config.h:
	cd $(BCFTOOLS_DIR); ./configure --disable-lzma
	touch $(BCFTOOLS_DIR)/config.h

$(SAMTOOLS_DIR):
	wget $(SAMTOOLS_URL) && tar -jxf $(SAMTOOLS_TAR) && rm -f $(SAMTOOLS_TAR)
	mv samtools-$(SAMTOOLS_VERSION) $(SAMTOOLS_DIR)

$(BCFTOOLS_DIR): 
	wget $(BCFTOOLS_URL) && tar -jxf $(BCFTOOLS_TAR) && rm -f $(BCFTOOLS_TAR)
	mv bcftools-$(BCFTOOLS_VERSION) $(BCFTOOLS_DIR)

$(UCSCTOOLS_DIR):
	wget $(UCSCTOOLS_URL) && tar -zxf $(UCSCTOOLS_TAR) && rm -f $(UCSCTOOLS_TAR)
	echo BINDIR=$(CURDIR)/bin > $(UCSCTOOLS_DIR)/kent/src/inc/localEnvironment.mk

$(UCSCTOOLS_LIB): $(UCSCTOOLS_DIR)
	cd $(UCSCTOOLS_DIR)/kent/src/lib && make

$(UCSCTOOLS_HTSLIB): $(UCSCTOOLS_DIR)
	cd $(UCSCTOOLS_DIR)/kent/src/htslib && make

$(BCFTOOLS_DIR)/plugins/%.c:
	ln -sf ../../gemBS_plugins/$(notdir $@) $(BCFTOOLS_DIR)/plugins/ 
	ln -sf ../../gemBS_plugins/$(basename $(notdir $@)).mk $(BCFTOOLS_DIR)/plugins/

$(FOLDER_BIN)/bedToBigBed: $(UCSCTOOLS_LIB) $(UCSCTOOLS_HTSLIB)
	cd $(UCSCTOOLS_DIR)/kent/src/utils/bedToBigBed && make

$(FOLDER_BIN)/wigToBigWig: $(UCSCTOOLS_LIB) $(UCSCTOOLS_HTSLIB)
	cd $(UCSCTOOLS_DIR)/kent/src/utils/wigToBigWig && make

ucsctools_clean:
	if [ -d $(UCSCTOOLS_DIR) ]; then $(MAKE) --directory=$(UCSCTOOLS_DIR)/kent/src/lib clean; $(MAKE) --directory=$(UCSCTOOLS_DIR)/kent/src/htslib clean; $(MAKE) --directory=$(UCSCTOOLS_DIR)/kent/src/utils/bedToBigBed clean; $(MAKE) --directory=$(UCSCTOOLS_DIR)/kent/src/utils/wigToBigWig clean; fi

clean: 	ucsctools_clean 
	@rm -f *~
	@rm -rf $(FOLDER_BIN) 
	if [ -f $(SAMTOOLS_DIR)/Makefile ]; then $(MAKE) --directory=$(SAMTOOLS_DIR) clean; fi
	if [ -f $(BCFTOOLS_DIR)/Makefile ]; then $(MAKE) --directory=$(BCFTOOLS_DIR) clean; fi
	if [ -f "gem3-mapper/Makefile.mk" ]; then $(MAKE) --directory=gem3-mapper clean; fi
	if [ -f "bs_call/src/Makefile.mk" ]; then $(MAKE) --directory=bs_call clean; fi

distclean: 	ucsctools_clean
	@rm -f *~
	@rm -rf $(FOLDER_BIN) 
	if [ -f $(SAMTOOLS_DIR)/Makefile ]; then cd $(SAMTOOLS_DIR); $(MAKE) clean; rm -f config.h config.log config.status config.mk; fi
	if [ -f $(BCFTOOLS_DIR)/Makefile ]; then cd $(BCFTOOLS_DIR); $(MAKE) clean; rm -f config.h config.log config.status config.mk; fi
	if [ -f "gem3-mapper/Makefile.mk" ]; then $(MAKE) --directory=gem3-mapper distclean; fi
	if [ -f "bs_call/src/Makefile.mk" ]; then $(MAKE) --directory=bs_call distclean; fi
