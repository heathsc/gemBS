GEMBS
=====

GEMBS is a bioinformatic pipeline designed for hightroughput analysis of DNA methylation from whole genome bisulfites sequencing data (WGBS). It implements GEM3, a high performance read aligner and BScall, a variant caller specifically for bisulfite sequencing data.


The manuscript describing the pipeline is available [here](https://www.biorxiv.org/content/early/2017/10/11/201988)

---------   
Licensing
---------

GEMBS is licensed under GPL.

--------
Download
--------

Use ``git clone --recursive`` to retrieve the complete source code including the code from external projects such as ``bs_call`` and ``gem3-mapper``.

    git clone --recursive https://github.com/heathsc/gemBS.git

------------
Installation
------------

1) Before starting the installation of GEMBS, please check if your system has the GSL library already installed.

    If your system does not have GSL library then you can download it from [GSL](https://www.gnu.org/software/gsl/) and follow the installation steps. 

    Once GSL is already available on your system then you can compile and install GEMBS.

2) Change GSL library paths. In order to compile bscall you must specify the GSL headers and library directories. 
   To do that, edit file ./tools/bs_call/Gsl.mk file with the proper paths. Just modify two lines starting with GSL_LIB and GSL_INC.

    ./tools/bs_call/Gsl.mk:  

        #1. MODIFY HERE THE GSL LIBRARY LOCATION. FOR Example: GSL_LIB = -L/path/to/gsl/lib
        GSL_LIB = -L/path/to/GSL/lib/
        #2. MODIFY HERE THE GSL HEADERS LOCATION. FOR Example: GSL_LIB = -L/path/to/gsl/include
        GSL_INC = -I/path/to/GSL/include/ 

3) Use python install command:

    ``python setup.py install --user``

-----------------------
Check your installation
-----------------------

For checking your installation follow the commands found at [README_TEST.md](test/README_TEST.md).


-------------
Documentation
-------------

Documentation can be found at [GEMBS](http://statgen.cnag.cat/GEMBS/)

----------
Changelog:
----------
    1.7 Fixed Fasta Contigs names with spaces in gem-index (gem3)
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
 
 GEMBS:
 * Marcos Fernandez-Callejo - marcos.fernandez@cnag.crg.eu
 
 gem mapper:
 * Santiago Marco-Sola - santiago.marco@cnag.crg.eu

 bisulfite caller and filtering:
 * Simon Heath - simon.heath@cnag.crg.eu


