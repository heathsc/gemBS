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
SAMTOOLS_VERSION=1.8
BCFTOOLS_VERSION=1.8
SAMTOOLS_DIR=samtools-$(SAMTOOLS_VERSION)
BCFTOOLS_DIR=bcftools-$(BCFTOOLS_VERSION)
SAMTOOLS=$(SAMTOOLS_DIR)/samtools
BCFTOOLS=$(BCFTOOLS_DIR)/bcftools
SAMTOOLS_TAR=$(SAMTOOLS_DIR).tar.bz2
BCFTOOLS_TAR=$(BCFTOOLS_DIR).tar.bz2
SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/$(SAMTOOLS_VERSION)/$(SAMTOOLS_TAR)
BCFTOOLS_URL=https://github.com/samtools/bcftools/releases/download/$(BCFTOOLS_VERSION)/$(BCFTOOLS_TAR)
INSTALL_PREFIX=/usr/local/gemBS

FOLDER_BIN=./bin/

all: release

static:	setup gem3-static $(SAMTOOLS) $(BCFTOOLS)
	$(MAKE) --directory=src static
	$(MAKE) --directory=gen static
	$(MAKE) --directory=filter_vcf static
	$(MAKE) --directory=cpgToWig static
	$(MAKE) --directory=sln static
	$(MAKE) --directory=bs_call static

release: setup	gem3 $(SAMTOOLS) $(BCFTOOLS)
	$(MAKE) --directory=src	
	$(MAKE) --directory=gen
	$(MAKE) --directory=filter_vcf
	$(MAKE) --directory=cpgToWig
	$(MAKE) --directory=sln
	$(MAKE) --directory=bs_call

debug: setup gem3-debug $(SAMTOOLS) $(BCFTOOLS)
	$(MAKE) --directory=src debug
	cd gen; ./configure
	$(MAKE) --directory=gen debug
	$(MAKE) --directory=filter_vcf debug
	$(MAKE) --directory=cpgToWig debug
	$(MAKE) --directory=sln debug
	$(MAKE) --directory=bs_call debug
	cd gem3-mapper; ./configure;

setup: gen/Makefile
	@mkdir -p $(FOLDER_BIN) 

gen/Makefile: gen/Makefile.in gen/configure
	cd gen; ./configure

gem3: gem3-mapper/Makefile.mk
	$(MAKE) --directory=gem3-mapper
	
gem3-static: gem3-mapper/Makefile.mk
	$(MAKE) --directory=gem3-mapper static
	
gem3-debug: gem3-mapper/Makefile.mk
	$(MAKE) --directory=gem3-mapper debug
	
gem3-mapper/Makefile.mk: gem3-mapper/Makefile.mk.in gem3-mapper/configure
	cd gem3-mapper; ./configure

$(SAMTOOLS_DIR)/config.mk: $(SAMTOOLS_DIR)
	cd $(SAMTOOLS_DIR); ./configure --prefix=$(INSTALL_PREFIX)

$(SAMTOOLS): $(SAMTOOLS_DIR)/config.h
	$(MAKE) --directory=$(SAMTOOLS_DIR) all all-htslib

$(BCFTOOLS_DIR)/config.h: $(BCFTOOLS_DIR)/plugins/mextr.c $(BCFTOOLS_DIR)/plugins/snpxtr.c
	cd $(BCFTOOLS_DIR); ./configure --prefix=$(INSTALL_PREFIX)

$(BCFTOOLS): $(BCFTOOLS_DIR)/config.h
	$(MAKE) --directory=$(BCFTOOLS_DIR) all

$(SAMTOOLS_DIR):
	wget $(SAMTOOLS_URL) && tar -jxf $(SAMTOOLS_TAR) && rm -f $(SAMTOOLS_TAR)
	
$(BCFTOOLS_DIR): 
	wget $(BCFTOOLS_URL) && tar -jxf $(BCFTOOLS_TAR) && rm -f $(BCFTOOLS_TAR)

$(BCFTOOLS_DIR)/plugins/%.c: $(BCFTOOLS_DIR)
	ln -sf ../../gemBS_plugins/$(notdir $@) $(BCFTOOLS_DIR)/plugins/ 
	ln -sf ../../gemBS_plugins/$(basename $(notdir $@)).mk $(BCFTOOLS_DIR)/plugins/

clean:
	@rm -f $(FOLDER_BIN)/* *~
	@rm -rf $(SAMTOOLS_DIR) $(BCFTOOLS_DIR)
	$(MAKE) --directory=gen clean
	$(MAKE) --directory=gem3-mapper clean
	$(MAKE) --directory=bs_call clean

