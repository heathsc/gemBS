gemBS
=====

gemBS is a wrapper to perform a set of different steps involved in the Bisulfite pipeline (Sequencing Mapping and Methylation Calling).

---------   
Licensing
---------

gemBS is licensed under GPL.

------------
Installation
------------

1) Before starting the installation of gemBS, you should check if your system has the GSL library already installed.

If your system does not have GSL library then you can download it from [GSL](https://www.gnu.org/software/gsl/) and follow the installation steps. 

Once GSL is already available on your system then you can compile and install gemBS.

2) Change GSL library paths. In order to compile bscall you must specify the GSL headers and library directories. 
   To do that, edit file ./tools/bs_call/Gsl.mk file with the proper paths. Just modify two lines starting with GSL_LIB and GSL_INC.

./tools/bs_call/Gsl.mk:

    #1. MODIFY HERE THE GSL LIBRARY LOCATION. FOR Example: GSL_LIB = -L/path/to/gsl/lib
    GSL_LIB = -L/path/to/GSL/lib/
    #2. MODIFY HERE THE GSL HEADERS LOCATION. FOR Example: GSL_LIB = -L/path/to/gsl/include
    GSL_INC = -I/path/to/GSL/include/ 

3) Use python install command:

    python setup.py install --user


-------------
Documentation
-------------

Documentation can be found at [gemBS] (http://statgen.cnag.cat/gemBS/)


----------
Changelog:
----------

    0.0 Initial Release  


-----------
Maintainers
-----------

Current maintainers:
 
 gemBS:
 * Marcos Fernandez-Callejo - marcos.fernandez@cnag.crg.eu
 
 gem mapper:
 * Santiago Marco-Sola - santiago.marco@cnag.crg.eu

 bisulfite caller and filtering:
 * Simon Heath - simon.heath@cnag.crg.eu


