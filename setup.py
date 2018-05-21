from distribute_setup import use_setuptools
use_setuptools()

import os
import sys
import glob

from distutils.core import setup
from setuptools import setup, Command
from setuptools.command.install import install as _install
from setuptools.command.build_py import build_py as _build_py

from distutils.command.clean import clean as _clean

import subprocess
import shutil

__VERSION_MAJOR = "2"
__VERSION_MINOR = "1"
__VERSION_SUBMINOR = "0"
__VERSION__ = "%s.%s.%s" % (__VERSION_MAJOR, __VERSION_MINOR,__VERSION_SUBMINOR)

def compile_gemBS_tools():
    process = subprocess.Popen(['make'], shell=True, cwd='tools')
    if process.wait() != 0:
        print >> sys.stderr, """
Error while compiling gemBS. That is very unfortunate.
A possible reason might be a missing dependency. Please take a look at the lines
before this one. You need the following programs and libraries installed to compile
the gemBS.
Programms needed:
    * make
    * gcc
Libraris needed:
    * python-dev (the python headers and include files)
    * libbz2-dev (for bz compression support)
On a Debian/Ubuntu system you should be able to get all needed dependencies with:
sudo apt-get install make gcc python-dev libbz2-dev
"""
        exit(1)


def clean_gemBS_tools():
    process = subprocess.Popen(['make clean'],shell=True,cwd='tools')
    if process.wait() != 0:
        print >> sys.stderr, """ Error Running cleaning. """
        exit(1)



def _install_bundle(install_dir):
    """Install compiled files to the GemBS installation Directory"""
    if install_dir is None:
        print "Unable to determine installation directory !"
        exit(1)
   
    gemBSbin_dir = "%s/%s" % (install_dir, "gemBSbinaries")
    if not os.path.exists(gemBSbin_dir):
        os.mkdir(gemBSbin_dir)

    # copy tools
    bins = [x for x in os.listdir("tools/bin")]
    for file in bins:
        if not file.endswith("gz"):
            print "Copy binary: %s to %s" % (file, gemBSbin_dir)
            result_file = "%s/%s" % (gemBSbin_dir, file)
            if os.path.exists(result_file):
                os.remove(result_file)
            shutil.copy("%s/%s" % ("tools/bin", file), gemBSbin_dir)
            os.chmod(result_file, 0755)

    # copy compiled bs_call tools
    bins = [x for x in os.listdir("tools/bs_call/bin")]
    for file in bins:
        print "Copy binary: %s to %s" % (file, gemBSbin_dir)
        result_file = "%s/%s" % (gemBSbin_dir, file)
        if os.path.exists(result_file):
            os.remove(result_file)
        shutil.copy("%s/%s" % ("tools/bs_call/bin", file), gemBSbin_dir)
        os.chmod(result_file, 0755)

    # copy compiled gem3 tools
    bins = [x for x in os.listdir("tools/gem3-mapper/bin")]
    for file in bins:
        print "Copy binary: %s to %s" % (file, gemBSbin_dir)
        result_file = "%s/%s" % (install_dir, file)
        if os.path.exists(result_file):
            os.remove(result_file)
        shutil.copy("%s/%s" % ("tools/gem3-mapper/bin", file), gemBSbin_dir)
        os.chmod(gemBSbin_dir, 0755)

    # copy samtools and bcftools
    bin_dir = "%s/%s" % (install_dir, "bin")
    lib_dir = "%s/%s" % (install_dir, "lib")
    plugins_dir = "%s/%s" % (install_dir, "libexec/bcftools")
    for dir in [bin_dir, lib_dir, plugins_dir]:
        if not os.path.exists(dir):
            os.makedirs(dir)
    if os.path.exists("tools/samtools/samtools"):
        print "Copy binary: samtools to {}".format(bin_dir)
        shutil.copy("tools/samtools/samtools", bin_dir)
    for htslib in glob.glob("tools/samtools/htslib*"):
        if os.path.isdir(htslib):
            for file in ["htsfile", "tabix", "bgzip"]:
                if os.path.exists("{}/{}".format(htslib,file)):
                    print "Copy binary: {} to {}".format(file, bin_dir)
                    shutil.copy("{}/{}".format(htslib,file), bin_dir)
            for file in ["libhts.a", "libhts.so"]:
                if os.path.exists("{}/{}".format(htslib,file)):
                    print "Copy library: {} to {}".format(file, lib_dir)
                    shutil.copy("{}/{}".format(htslib,file), lib_dir)
                    
            
    if os.path.exists("tools/bcftools/bcftools"):
        print "Copy binary: bcftools to {}".format(bin_dir)
        shutil.copy("tools/bcftools/bcftools", bin_dir)
    plugins = [x for x in glob.glob("tools/bcftools/plugins/*.so")]
    for file in plugins:
        print "Copy plugin: {} to {}".format(file, plugins_dir)
        shutil.copy(file, plugins_dir)

# hack the setup tools installation
class install(_install):

    def run(self):
        _install.run(self)
        
        # find target folder
        install_dir = None
        for file in self.get_outputs():
           if file.endswith("/gemBS/__init__.py"):
                install_dir = "%s" % os.path.split(file)[0]
                break

        _install_bundle(install_dir)
        
 
# hack the setup tools building       
class build_py(_build_py):
    
    def run(self):
        compile_gemBS_tools()
        parent_dir = os.path.split(os.path.abspath(__file__))[0]
        target_dir = "%s/%s" % (parent_dir, "gemBS")
        _install_bundle(target_dir)
        _build_py.run(self)
        
# hack the setup tools cleaning
class clean(_clean):

    def run(self):
        _clean.run(self)
        clean_gemBS_tools()


_commands = {'install': install,'build_py': build_py,'clean':clean}


setup(cmdclass=_commands,
      name='gemBS',
      version=__VERSION__,
      description='Python application to perform the different steps involved in the Bisulphite Pipeline.',
      author='Marcos Fernandez-Callejo, Santiago Marco-Sola, Simon Heath',
      author_email='marcos.fernandez@cnag.crg.eu',
      url='http://statgen.cnag.cat/gemBS/',
      packages=['gemBS'],
      package_data={"": ["%s/%s" % ("gemBS/gemBSbinaries", x) for x in ["readNameClean",
                                                                      "filter_vcf",
                                                                      "cpgToWig",
                                                                      "align_stats",
                                                                      "gem-constructor",
                                                                      "gem-indexer",
                                                                      "gem-mapper",
                                                                      "gem-retriever",
                                                                      "bs_call",
                                                                      "dbSNP_idx"
                                                                     ]]},
      entry_points = {
        'console_scripts': ['gemBS=gemBS.commands:gemBS_main'],
      }
     )

