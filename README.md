GEMBS
=======

GEMBS is a bioinformatic pipeline designed for hightroughput analysis of DNA methylation from whole genome bisulfites sequencing data (WGBS). It implements GEM3, a high performance read aligner and BScall, a variant caller specifically for bisulfite sequencing data.


The manuscript describing the pipeline is available [here](http://www.webmanuscript.com)

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


-------------
Documentation
-------------

Documentation can be found at [GEMBS](http://statgen.cnag.cat/gemBS/)


----------
Changelog:
----------

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


