bscall
======

Variant Caller for Bisulfite Sequencing Data.


------------
Installation
------------

If your system does not have the hts library (htslib - used by samtools and
bcftools) then you should install it before proceeding.

Configure:

The compilation should be configured by typing:

	./configure
	
If htslib is not in a standard location on your system then the location
should be specified using the --with-htslib option to configure. For
example, if the installation prefix for htslib is /opt/local
(so the libraries can be found in /opt/local/lib and the include
directory htslib in /opt/local/include) then the configuration command
line should be:

	./configure --with-htslib=/opt/local

Compile:

    make all

Install:

If the compilation process has been successfully completed then a
binary file should be found at bin directory. Just copy it to a
directory included in your $PATH.

Copy binary:

    cp bin/bscall /usr/bin

--------------
Running bscall
--------------

Run bscall from a BAM file to get a bcf output:

    bs_call -r my_reference.fasta -p -L5 -n my_sample_name -o mysample.bcf mysample.bam

The parameters configured for this example are -p (Paired End Data) and -L5 (5 bases to trim from left of read pair).

---------
Changelog
---------
    2.1.7 Fix incorrect range check
    2.1.6 Add option for read end specific trimming.  Add debug target to Makefile.  
          Add range checking to release target. Make changes for htslib 1.11.
    2.1.5 Move to new dbSNP format, allowing loading of individual chromosomes
    2.1.4 Fix bug when using the -k flag, and fix a rare hang at end of input
    2.1.3 If output format is not explicilty specified, an attempt to guess the format from the filename is made (i.e.,
	       if the output file ends in .bcf or .vcf.gz the appropriate file type will be selected).
	 2.1.3 If not output name is set (so output is to standard out) and standard out is a terminal, output format is set 
	       to uncompressed VCF, irrespective of the specified output format.
    2.1.3 Removed calculation of GOF (goodness of fit) filter as it is expensive, non-standard and unreliable
    2.1.2 Improved the support for CRAM input
    2.1.1 Switched to htslib 1.10
          Added --benchmark-mode where date and version numbers are not written to output files
    2.1.0 Reorganized and cleaned up distribution.  
          Switched to using htslib for input and output, so can now read from SAM or BAM and
          write to VCF or BCF natively.
          Reduced *a lot* the memory usage by (a) only reading the reference for the contig 
          currently being processed and (b) reducing the amount of time most alignment
          information is kept.  A typical human WGBS sample can no be processed calling
          all chromosomes in parallel on a single computer using < 10GB RAM (this will
          depend on the coverage).
          Changed the threading model.  Additional threads are split between calculation,
          input and ouput.
          Removed most of the unused GEMTools library (as we no longer use it for parsing the 
          SAM input, keeping just the core library which can now be found in gt/
          Removed legacy and unused options.
    2.0.3 Add distclean target to makefile
    2.0.3 Remove compile warnings from GEMTools
    2.0.3 Fix bug with handling duplicate reads  
    2.0.2 Document configuration process.
    2.0.1 Fix argument -k about discarded reads that do not form proper pairs.
    2.0.1 Fix Single End Memory Leak.
    2.0.1 Use of dbSNP to evaluate SNP Calling.
    2.0.1 Output JSON stats to measure SNP and Methylation calls.
    2.0.1 Index set of Bed files for dbSNP.

---------
Licensing
---------

bscall is licensed under GPL. See LICENSE for more information.

------
Author
------

Simon Heath at (CNAG/CRG) Centre Nacional d’Anàlisi Genòmica / Centre de Regulació Genòmica.
simon.heath@gmail.com
