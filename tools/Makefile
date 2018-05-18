#============================================================================
# PROJECT: gemBS
# FILE: Makefile
# DATE: 01/01/2015
# AUTHOR(S): Marcos Fernandez Callejo <mfernandez@cnag.crg.eu>
# DESCRIPTION: Top level makefile
#============================================================================

# Definitions
ROOT_PATH=$(CURDIR)

FOLDER_BIN=./bin/

all: release

static:	setup gem3-static
	$(MAKE) --directory=src static
	$(MAKE) --directory=gen static
	$(MAKE) --directory=filter_vcf static
	$(MAKE) --directory=cpgToWig static
	$(MAKE) --directory=sln static
	$(MAKE) --directory=bs_call static

release: setup	gem3
	$(MAKE) --directory=src	
	$(MAKE) --directory=gen
	$(MAKE) --directory=filter_vcf
	$(MAKE) --directory=cpgToWig
	$(MAKE) --directory=sln
	$(MAKE) --directory=bs_call

debug: setup gem3-debug
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
	
clean:
	@rm -f $(FOLDER_BIN)/* *~
	$(MAKE) --directory=gen clean
	$(MAKE) --directory=gem3-mapper clean
	$(MAKE) --directory=bs_call clean

