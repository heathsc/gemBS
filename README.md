gemBS
=====

gemBS is a high performance bioinformatic pipeline designed for highthroughput analysis
of DNA methylation data from whole genome bisulfites sequencing data
(WGBS). It combines GEM3, a high performance read aligner and
bs_call, a high performance variant and methyation caller, into a streamlined and efficient pipeline for
bisulfite sueqnce analysis.

The manuscript describing the pipeline is available [here](https://www.biorxiv.org/content/early/2017/10/11/201988)

---------   
Licensing
---------

gemBS is licensed under GPL.

--------
Download
--------

Use ``git clone --recursive`` to retrieve the complete source code including the code from external projects such as ``bs_call`` and ``gem3-mapper``.

    git clone --recursive https://github.com/heathsc/gemBS.git

------------
Installation
------------

1) Before starting the installation of gemBS, you will need to install
or check the installation of several packages.

  a) gcc with development libraries
  b) python3, pip3, matplotlib, multiprocess
  c) zlib, gsl, libncurses
  
If you are working on a clean (fairly recent) Ubuntu installation, you
can install everything required with the followiwg commands:

    sudo apt-get update
    sudo apt-get install -y python3 build-essential git python3-pip wget pigz
    sudo apt-get install -y zlib1g-dev libbz2-dev gsl-bin libgsl0-dev
    sudo apt-get install -y libncurses5-dev liblzma-dev libssl-dev libcurl4-openssl-dev
    pip3 install matplotlib multiprocess

2) Download the gemBS distribution if you haven't already done so:

    ``git clone --recursive https://github.com/heathsc/gemBS.git``

3) Use python install command:

To install to the standard system location (i.e., so that all users
can use gemBS):

    ``python3 setup.py install``

To install to the user's home directory:

    ``python3 setup.py install --user``

Note that gemBS requires that GSL (GNU scientific library) is
installed prior to the installation of gemBS.  If GSL has been
installed on your system to the standard system location then the
above procedure should work without modification.  If, however, GSL
has been installed to a non-standard location then then --gsl-path
option to install can be used.  Fo example, if the installation prefix
for GSL is /opt/local (so the libraries are in /opt/local/lib and the
include directory gsl is in /opt/local/include) then the following
install command should be used:

    ``python3 setup.py install --gsl-path=/opt/local``

or

    ``python3 setup.py install --gsl-path=/opt/local --user``

-----------------------
Check your installation
-----------------------

For checking your installation follow this
[worked example](http://statgen.cnag.cat/gemBS/UserGuide/_build/html/example.html).


-------------
Documentation
-------------

Documentation can be found at
[gemBS](http://statgen.cnag.cat/gemBS/)

----------
Changelog:
----------
    3.4.2 Add CRAM support (via make_cram option in configuration file)
    3.4.1 Add benchmark-mode that does not write date or program version numbers into SAM/BAM or VCF/BCF files
	       Switch to samtools, bcftools and htslib v1.10
    3.4.0 Move to new bs_call version (2.1.0) which is more efficient
          in memory use and can read BAMs and write BCFs natively.
          The new bs_call requires a faidx indexed reference, so gemBS
          no creates this during indexing.
	 3.4.0 Add switches to give more control to threads and memory
          usage in mapping and calling stages
    3.3.3 Remove legacy pathway for config files with no header line (fix issue 'error in gemBS index #65)
    3.3.2 Fix error where header line for wig files could be omitted
    3.3.2 Fix generation of non_cpg files
    3.3.1 Fix Attribute error bug due to not checking if conversion is a list
    3.3.0 Make new release for IHEC
    3.3.0 Switch conversion default in IHEC_standard configuration to 0.01,0.05 rather than auto, which can give odd results if conversion controls not present or not working correctly
    3.3.0 Fix bug where conversion parameters could be ignored
    3.2.13 Fix formatting bug in mextr with multiple samples (not triggered in normal gemBS use)
    3.2.12 Ensure that conversion statistics are correctly calculated for non-stranded or reverse conversion protocols
    3.2.11 Introduce reverse_conversion option for mapping where read 1 is G2A converted and read 2 is C2T converted
    3.2.10 Correct regex patch for single end reads
    3.2.9 Update Singularity and Dockerfile recipes to allow kemp utils to be built correctly
    3.2.9 Make setup.py and gemBS/commands.py read the version information from gemBS/version.py, so ensuring consistency
    3.2.9 Fix bug added in last version where options in config file were not being taken into account
    3.2.8 Fix mis specification errors in long options for mextr. 
    3.2.8 Fix bug where mextr (methyl extract plugin for bcftools) would crash if cpg output  option was not set.
    3.2.7 Apply patches for bugs in handling single end reads (suggested by I. Moghul)
    3.2.7 Changed regex for filenames to make it more general (suggested by I. Moghul)
    3.2.7 Fixed bug (reported by chhylp123) where zero arguments to some options were being ignored
    3.2.6 Cleaned up compilation and cleaning of gemBS tools
    3.2.6 Fixed python error if either the over conversion reference sequence was not defined
    3.2.6 Added check in bs_call that conversion parameters are valid (between 0 and 1)
    3.2.6 Perform more stringent sanity checking on conversion vaalues when autocomputed by gemBS
    3.2.6 Use --diasble-lzma configuration flag for samtools and bcftools as we don't need it and it removes an unneccesary dependency
    3.2.6 Add install options --disable-cuda (on by default) and --enable-cuda that affect GEM3 comppilation
    3.2.6 Bug fix with incorrect handling of duplicate reads
    3.2.5 Minor bug fix - correct error with non-paired end non-bisulfite reads
    3.2.4 Modify the bisulfite processing in gem-mapper to be more efficient (in particular for the non-stranded option)
    3.2.4 Modify gemBS to use the new conversion options for gem-mapper
    3.2.4 Switch gem-mapper to use option --underconversion-sequence and --overconversion-sequence rather than --underconversion_sequence to be consistent with other options
    3.2.3 Fixed bug if conversion parameters were not set
    3.2.2 Rework non-stranded mode so that both possible conversions are tried and the results merged
    3.2.2 Fix bug where non-stranded flag was not being passed to mapper in paired end mode
    3.2.1 Move warning message from bscall from stdout to stderr
    3.2.1 Switch Singularity build to use Ubuntu 16.04 rather than 18.04 to allow the image to work in CentOS 6 (Docker build changed as well to keep the two in sync)
    3.2.1 Fix undeclared variable bugs and missing --ignore-deps option in merge-bcfs
    3.2.1 Add default for dbSNP_index if dbSNP_files is set
    3.2.1 Add gsl-path install option
    3.2.0 Make new release
    3.1.0 Make installation process more modular.  Allow for sub-installs
    3.1.0 Add support for reading config from ${index_dir}/gemBS.json if it exists
    3.1.0 Add --reference-bias option to mextr and gemBS extract
    3.1.0 Add support for non-bisulfite mapping of individual datasets
    3.1.0 Allow white space in variable values
    3.1.0 Allow fallback to gzip if pigz not present
    3.1.0 Add --dry-run, --json, --ignore-db and --ignore-dep to extract command
    3.1.0 Add --ignore-dep option to call and merge-bcfs commands
    3.1.0 Add SNP extraction function to extract command
    3.0 Make v3.0 release
    3.0 Merge with master branch.
    3.0 Bump samtools sort memory limit to 2G
    3.0 Add extra_references option for reference generation
    3.0 Allow input files to mapping to be shell commands
    3.0 Add links to documentation
    3.0 Upload new yeast example and add documentation
    3.0 Add --dir option to gemBS
    3.0 Add --ignore-db options for --dry-run / --json
    3.0 Add --json output option for dry runs
    3.0 Update help text to match new functions
    3.0 Introduce standard analysis configurations stored within distribution
    3.0 Switch gem3-mapper distribution to gembs branch on official gem3-mapper repo
    3.0 Removal of incomplete files and roll back of db in the event of pipeline failure
    3.0 Automatic removal of individual BAMs and BCFs after successful merging
    3.0 Prevent pipelines hanging in event of failure
    3.0 Generate ENCODE bed and bigbed files
    3.0 Switch to python 3
    3.0 Switch to mextr for BCF filtering
    3.0 Include fetch and build of samtools / bcftools during build process
    3.0 Add dry-run capability to map and call commands
    3.0 Introduce contig pools to automatically group small contigs
    3.0 Automatic generation of contig.size files from index build
    3.0 Allow use of in memory sqlite3 db as an option
    3.0 Allow multiple instances of gemBS (possible on different hosts) to work 
        simultaneously on the same analysis
    3.0 Reduce and simply commands
    3.0 Add Dockerfile
    3.0 Add multi-threading and multi-processing options for most commands
    3.0 Use sqlite3 to track progress of analyses, file paths etc.
    3.0 Added more flexible configuration options (new csv format + new configuration file)
    3.0 Remove test dataset from distribution (distribute from web site)
    2.1.0 Ensure commands run during pipeline come from installation
    2.1.0 Added Singularity build recipe
    2.1.0 Add new command gemBS direct-mapping
    2.1.0 Fixed Makefile clean in tools
    2.0.2 Fixed bug related with the percentage of High Quality Variant in Variants summary report.
    2.0.2 Check temporary directory existence.
    2.0.2 Fixed QualityNonRefCpg sample name in png image.
    2.0.2 Fixed mapper issues related with aligning performace.
    2.0.2 Fixed arguments for Under/Over Conversion sequence name in gem3-mapper
    2.0.1 On bscall repository, fixed argument -k about discarded reads that do not form proper pairs.
    2.0 Check tmp folder before starting mapping process.
    2.0 Added Left and Right Trimming optional arguments to gemBS bscall.
    2.0 Added GC Coverage correlation value to BS Call Stats Summary.
    2.0 Fixed error when reporting complete path to not found bam files.
    2.0 Fixed iteration over sampleBams dictionary in MergeAll method.
    2.0 Updated: Avoid redo indexing when merging just one file.
    2.0 Changed conversion formula.
    2.0 Added parameter for dbSNP.
    2.0 Added threads to bscall.
    2.0 Removed CpGs reports. Already done from bscall report.
    2.0 Fixed bs_call makefile for the gcc to be used.
    2.0 New bscall version. Generates JSON report.
    2.0 Removed gemBS options snp-stats,cpg-report,cpg-stats.
    2.0 Added summary report from the bs_call json stats
    2.0 New BSCall Report. From bscall son file generates three types of reports:
        Mapping and Coverage Report
        Bs-Genotypes Calls Report
        Methylation Statistics report
    1.7 Added non stranded read conversion parameter
    1.7 Fixed SE crash when estimating overlapped bases.
    1.7 Fixed gem-index (gem3) to follow fastq and SAM specifications. 
        Modified gem3-mapper repository external module.
        New external module https://github.com/heathsc/gem3-mapper.git
    1.7 Fixed threads parameter to samtools merge
    1.7 Fixed threads parameter to gem-mapper
    1.7 Removed Indels Field on Variants Report.
    1.7 Added Overlapping Bases at Mapping Report
    1.7 Modified Base Counts Overall, removed Base Counts general and Base Counts Overall
    1.7 New Dinucleotide CpGs Report
        New table dinucleotide stats
        New plots for Informative Reads and CpGs
        Methylation levels plots for different types of CpGs
        Summary Table
    1.7 New Readme file to inform about report test
    1.7 New basic statis table for Variants Report
    1.7 Removed parameter -r (reference length) parameter for mapping reports command (gemBS bsMap).
    1.6 New CpGs Density plot, include box plos, bar plot and fitting curve
    1.6 Change name at CpG report:
        "Heterozygous" for "Alternative CX"
        "De Novo CpGs Methylation Status" for "Non Reference CpGs"
        "CpGs with SNP" for "SNPs (CX) at Reference CpGs"
    1.6 CpGs Report Simplified to Q>20
    1.6 BigWig Default parameters for filtering CpG per a given quality and a total number of supported informative reads   
    1.5 Initial Release  


----------
Developers
----------
 
 gemBS:
 * Marcos Fernandez-Callejo - marcos.fernandez@cnag.crg.eu
 * Simon Heath - simon.heath@gmail.com
 
 gem mapper:
 * Santiago Marco-Sola - santiagomsola@gmail.com

 bisulfite caller and filtering:
 * Simon Heath - simon.heath@gmail.com


