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

    ``python3 setup.py install --user``

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


