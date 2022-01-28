#============================================================================
# PROJECT: bscall
# FILE: Makefile
# DATE: 02/05/2017
# AUTHOR(S): Marcos Fernandez Callejo <mfernandez@cnag.crg.eu>
# DESCRIPTION: Top level makefile
#============================================================================

# Definitions
all: release

release: 	
	$(MAKE) --directory=gt
	$(MAKE) --directory=src

debug:
	$(MAKE) --directory=gt debug
	$(MAKE) --directory=src debug

static:
	$(MAKE) --directory=gt static
	$(MAKE) --directory=src static

clean:
	$(MAKE) --directory=gt clean
	$(MAKE) --directory=src clean

distclean:
	$(MAKE) --directory=gt distclean
	$(MAKE) --directory=src distclean
	rm -f config.status config.log
