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
SAMTOOLS_VERSION=1.10
BCFTOOLS_VERSION=1.10.2
HTSLIB_VERSION=1.10.2
SAMTOOLS_DIR=samtools
BCFTOOLS_DIR=bcftools
HTSLIB_DIR=htslib
SAMTOOLS=$(SAMTOOLS_DIR)/samtools
BCFTOOLS=$(BCFTOOLS_DIR)/bcftools
SAMTOOLS_TAR=samtools.tar.gz
BCFTOOLS_TAR=bcftools.tar.gz
HTSLIB_TAR=htslib.tar.gz
SAMTOOLS_URL=https://github.com/samtools/samtools/archive/$(SAMTOOLS_VERSION).tar.gz
BCFTOOLS_URL=https://github.com/samtools/bcftools/archive/$(BCFTOOLS_VERSION).tar.gz
HTSLIB_URL=https://github.com/samtools/htslib/archive/$(HTSLIB_VERSION).tar.gz

MACHTYPE:=$(shell uname -m)
ifneq (,$(findstring -,$(MACHTYPE)))
    MACHTYPE:=$(shell uname -m)
endif
		  
FOLDER_BIN=bin

all: setup gem3 _samtools _bcftools _bs_call _utils

all_static: setup gem3-static _samtools _bcftools _bs_call _utils

all_debug: setup gem3-debug _samtools _bcftools _bs_call _utils

_samtools: $(HTSLIB_DIR) $(SAMTOOLS_DIR) $(SAMTOOLS_DIR)/config.mk 
	$(MAKE) --directory=$(SAMTOOLS_DIR) all all-htslib

_bcftools: $(HTSLIB_DIR) $(BCFTOOLS_DIR) $(BCFTOOLS_DIR)/config.h
	$(MAKE) $(FOLDER_BIN)/bcftools

_bs_call: bs_call/src/Makefile.mk
	$(MAKE) --directory=bs_call

_utils: utils/Makefile 
	$(MAKE) --directory=utils

utils/Makefile: utils/Makefile.in utils/configure _samtools
	cd utils; ./configure --with-htslib=../${HTSLIB_DIR}

setup:
	@mkdir -p $(FOLDER_BIN) 

bs_call/src/Makefile.mk: bs_call/src/Makefile.mk.in bs_call/gt/Makefile.mk.in bs_call/configure _samtools
	cd bs_call; ./configure ${BS_CALL_CONFIG} --with-htslib=../../${HTSLIB_DIR}

gem3: gem3-mapper/Makefile.mk
	$(MAKE) --directory=gem3-mapper
	
gem3-static: gem3-mapper/Makefile.mk
	$(MAKE) --directory=gem3-mapper static
	
gem3-debug: gem3-mapper/Makefile.mk
	$(MAKE) --directory=gem3-mapper debug
	
gem3-mapper/Makefile.mk: gem3-mapper/Makefile.mk.in gem3-mapper/configure
	cd gem3-mapper; ./configure ${GEM3_CONFIG}

$(SAMTOOLS_DIR)/config.mk:
	cd $(SAMTOOLS_DIR); ./configure # --disable-lzma

$(FOLDER_BIN)/bcftools: $(BCFTOOLS_DIR) 
	$(MAKE) --directory=$(BCFTOOLS_DIR) all

$(BCFTOOLS_DIR)/config.h:
	cd $(BCFTOOLS_DIR); ./configure # --disable-lzma
	touch $(BCFTOOLS_DIR)/config.h

$(HTSLIB_DIR):
	wget -O $(HTSLIB_TAR) $(HTSLIB_URL) && tar -zxf $(HTSLIB_TAR) && rm -f $(HTSLIB_TAR)
	mv htslib-$(HTSLIB_VERSION) $(HTSLIB_DIR)
	(cd $(HTSLIB_DIR) && cat ../patches/*.patch | patch && autoheader && autoconf && ./configure)

$(SAMTOOLS_DIR):
	wget -O $(SAMTOOLS_TAR) $(SAMTOOLS_URL) && tar -zxf $(SAMTOOLS_TAR) && rm -f $(SAMTOOLS_TAR)
	mv samtools-$(SAMTOOLS_VERSION) $(SAMTOOLS_DIR)
	(cd $(SAMTOOLS_DIR) && autoreconf)

$(BCFTOOLS_DIR): 
	wget -O $(BCFTOOLS_TAR) $(BCFTOOLS_URL) && tar -zxf $(BCFTOOLS_TAR) && rm -f $(BCFTOOLS_TAR)
	mv bcftools-$(BCFTOOLS_VERSION) $(BCFTOOLS_DIR)
	(cd $(BCFTOOLS_DIR) && autoheader && autoconf)

clean:
	@rm -f *~
	@rm -rf $(FOLDER_BIN) 
	if [ -f $(SAMTOOLS_DIR)/Makefile ]; then $(MAKE) --directory=$(SAMTOOLS_DIR) clean; fi
	if [ -f $(BCFTOOLS_DIR)/Makefile ]; then $(MAKE) --directory=$(BCFTOOLS_DIR) clean; fi
	if [ -f "gem3-mapper/Makefile.mk" ]; then $(MAKE) --directory=gem3-mapper clean; fi
	if [ -f "bs_call/src/Makefile.mk" ]; then $(MAKE) --directory=bs_call clean; fi
	if [ -f "utils/Makefile" ]; then $(MAKE) --directory=utils clean; fi

distclean:
	@rm -f *~
	@rm -rf $(FOLDER_BIN) 
	if [ -f $(SAMTOOLS_DIR)/Makefile ]; then cd $(SAMTOOLS_DIR); $(MAKE) clean; rm -f config.h config.log config.status config.mk; fi
	if [ -f $(BCFTOOLS_DIR)/Makefile ]; then cd $(BCFTOOLS_DIR); $(MAKE) clean; rm -f config.h config.log config.status config.mk; fi
	if [ -f $(HTSLIB_DIR)/Makefile ]; then cd $(HTSLIB_DIR); $(MAKE) clean; rm -f config.h config.log config.status config.mk; fi
	if [ -f "gem3-mapper/Makefile.mk" ]; then $(MAKE) --directory=gem3-mapper distclean; fi
	if [ -f "bs_call/src/Makefile.mk" ]; then $(MAKE) --directory=bs_call distclean; fi
	if [ -f "utils/Makefile" ]; then $(MAKE) --directory=utils distclean; fi
