import os
import sys
import glob

from distutils.core import setup
from distutils.sysconfig import get_config_vars
from setuptools import setup, Command
from setuptools.command.install import install as _install
# from setuptools.command.build_py import build_py as _build_py
from setuptools.command.install_scripts import install_scripts as _install_scripts
from setuptools.command.install_lib import install_lib as _install_lib
from distutils.command.clean import clean as _clean

from site import USER_BASE
from site import USER_SITE

import subprocess
import shutil

pwd = os.path.abspath(os.path.dirname(__file__))
exec(open(os.path.join(pwd,'gemBS','version.py')).read())

def compile_gemBS_tools(options, gsl_path, enable_cuda, disable_cuda):
    make_com = 'make ' + ' '.join(options)
    if gsl_path:
        make_com = "BS_CALL_CONFIG=\'--with-gsl={}\' ".format(gsl_path) + make_com
    if disable_cuda:
        make_com = "GEM3_CONFIG=\'--disable-cuda\' " + make_com
    elif enable_cuda:
        make_com = "GEM3_CONFIG=\'--enable-cuda\' " + make_com
    process = subprocess.Popen(make_com, shell=True, cwd='tools')
    if process.wait() != 0:
        print ("""
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
""", file=sys.stderr)
        exit(1)


def clean_gemBS_tools():
    process = subprocess.Popen(['make distclean'],shell=True,cwd='tools')
    if process.wait() != 0:
        print (""" Error Running cleaning. """,file=sys.stderr)
        exit(1)



def _install_bundle(install_dir, inst):
    """Install compiled files to the GemBS installation Directory"""
    if install_dir is None:
        print("Unable to determine installation directory !")
        exit(1)
   
    gemBSbin_dir = os.path.join(install_dir, "gemBSbinaries")
    if not os.path.exists(gemBSbin_dir):
        os.mkdir(gemBSbin_dir)

    # copy tools/bin
    bins = ['gemBS_cat', 'readNameClean', 'md5_fasta']
    if not (inst.minimal or inst.no_kent):
        bins.extend(['wigToBigWig', 'bedToBigBed'])
    for file in bins:
        f = os.path.join('tools/bin', file)
        if os.path.exists(f):
            # print ("Copy binary: %s to %s" % (file, gemBSbin_dir))
            result_file = os.path.join(gemBSbin_dir, file)
            if os.path.exists(result_file):
                os.remove(result_file)
            shutil.copy(f, gemBSbin_dir)
            os.chmod(result_file, 0o755)

    if not (inst.minimal or inst.no_bscall):
        # copy compiled bs_call tools
        bins = [x for x in os.listdir("tools/bs_call/bin")]
        for file in bins:
            # print ("Copy binary: %s to %s" % (file, gemBSbin_dir))
            result_file = os.path.join(gemBSbin_dir, file)
            if os.path.exists(result_file):
                os.remove(result_file)
            shutil.copy(os.path.join("tools/bs_call/bin", file), gemBSbin_dir)
            os.chmod(result_file, 0o755)

    if not (inst.minimal or inst.no_gem3):
        # copy compiled gem3 tools
        bins = [x for x in os.listdir("tools/gem3-mapper/bin")]
        for file in bins:
            # print ("Copy binary: %s to %s" % (file, gemBSbin_dir))
            result_file = os.path.join(install_dir, file)
            if os.path.exists(result_file):
                os.remove(result_file)
            shutil.copy(os.path.join("tools/gem3-mapper/bin", file), gemBSbin_dir)
            os.chmod(gemBSbin_dir, 0o755)


    # copy samtools, bcftools and config files
    bin_dir = os.path.join(install_dir, "bin")
    lib_dir = os.path.join(install_dir, "lib")
    plugins_dir = os.path.join(install_dir, "libexec", "bcftools")
    etc_dir = os.path.join(install_dir, "etc")
    config_dir = os.path.join(etc_dir, "gemBS_configs")
    for dir in [bin_dir, lib_dir, plugins_dir, config_dir]:
        if not os.path.exists(dir):
            os.makedirs(dir)
    if not (inst.minimal or inst.no_samtools):
        if os.path.exists("tools/samtools/samtools"):
            # print ("Copy binary: samtools to {}".format(bin_dir))
            shutil.copy("tools/samtools/samtools", bin_dir)
            os.chmod(os.path.join(bin_dir, "samtools"), 0o755)
        for htslib in glob.glob("tools/samtools/htslib*"):
            if os.path.isdir(htslib):
                for file in ["htsfile", "tabix", "bgzip"]:
                    if os.path.exists(os.path.join(htslib,file)):
                        # print ("Copy binary: {} to {}".format(file, bin_dir))
                        shutil.copy(os.path.join(htslib,file), bin_dir)
                        os.chmod(os.path.join(bin_dir, file), 0o755)
                for file in ["libhts.a", "libhts.so", "libhts.dylib"]:
                    if os.path.exists(os.path.join(htslib,file)):
                        # print ("Copy library: {} to {}".format(file, lib_dir))
                        shutil.copy(os.path.join(htslib,file), lib_dir)
                        os.chmod(os.path.join(lib_dir, file), 0o755)
                                
    if os.path.exists("tools/bcftools/bcftools"):
        # print ("Copy binary: bcftools to {}".format(bin_dir))
        shutil.copy("tools/bcftools/bcftools", bin_dir)
        os.chmod(os.path.join(bin_dir, "bcftools"), 0o755)
    plugins = [x for x in glob.glob("tools/bcftools/plugins/*.so")]
    for file in plugins:
        # print ("Copy plugin: {} to {}".format(file, plugins_dir))
        shutil.copy(file, plugins_dir)
        os.chmod(os.path.join(plugins_dir,os.path.basename(file)), 0o755)

    files = [x for x in os.listdir("gemBS/etc")]
    for file in files:
        if os.path.isfile(os.path.join("gemBS/etc",file)):
            # print ("Copy {} to {}".format(file, etc_dir))
            shutil.copy(os.path.join("gemBS/etc",file), etc_dir)
            os.chmod(os.path.join(etc_dir, file), 0o644)
    
    files = [x for x in os.listdir("gemBS/etc/gemBS_configs")]
    for file in files:
        # print ("Copy {} to {}".format(file, config_dir))
        shutil.copy(os.path.join("gemBS/etc/gemBS_configs",file), config_dir)
        os.chmod(os.path.join(config_dir, file), 0o644)
        
# hack the setup tools installation
class install(_install):
    _install.user_options.extend([
        ('no-samtools', None, "Do not install samtools"),
        ('no-kent', None, "Do not install kent tools (bedToBigBed, wigToBigWig)"),
        ('no-gem3', None, "Do not install gem3 mapper"),
        ('no-bscall', None, "Do not install bscall"),
        ('minimal', None,
         "Perform minimal install (equivalent to --no-samtools --no-kent --no-gem3 --no-bscall)"),
        ('disable-cuda', None, "Do not build GPU support for GEM3 (default)"),
        ('enable-cuda', None, "Try to build GPU support for GEM3"),
        ('gsl-path=', None, "Installation path of GSL library")
    ])
    _install.boolean_options.extend(['no-samtools','no-kent','no-gem3','no-bscall','minimal'])

    def initialize_options(self):
        self.minimal = False
        self.no_samtools = False
        self.no_kent = False
        self.no_gem3 = False
        self.no_bscall = False
        self.disable_cuda = False
        self.enable_cuda = False
        self.gsl_path = None
        _install.initialize_options(self)
        
    def run(self):
        options=['setup', '_bcftools', '_utils']
        if not self.minimal:
            if not self.no_samtools:
                options.append('_samtools')
            if not self.no_kent:
                options.append('kent')
            if not self.no_gem3:
                options.append('gem3')
            if not self.no_bscall:
                options.append('_bs_call')

        if not self.enable_cuda:
            self.diable_cuda = True
            
        compile_gemBS_tools(options, self.gsl_path, self.enable_cuda, self.disable_cuda)
        _install.run(self)

        # find target folder
        install_dir = os.path.join(gemBS_install_dir, "gemBS")
        _install_bundle(install_dir, self)
        
 
class install_lib(_install_lib):
    def run(self):
        # Store install_dir for future use
        global gemBS_install_dir
        gemBS_install_dir = self.install_dir
        _install_lib.run(self)
    
    def get_install_dir(self):
        return self.install_dir
    
class install_scripts(_install_scripts):
    
    def write_script(self, script_name, contents, mode="t", *ignored):
        i = contents.find('__requires__')
        if i >= 0:
            j = contents.rfind('\n', 0, i)
            if j >= 0:
                contents = contents[:j+1] + "import sys\nsys.path.insert(0,'{}')\n".format(gemBS_install_dir) + contents[j+1:]
        _install_scripts.write_script(self, script_name, contents, mode, *ignored)
            
# hack the setup tools cleaning
class clean(_clean):

    def run(self):
        _clean.run(self)
        clean_gemBS_tools()


#_commands = {'install': install,'install_lib': install_lib, 'install_scripts': install_scripts, 'build_py': build_py,'clean':clean}
_commands = {'install': install,'install_lib': install_lib, 'install_scripts': install_scripts,'clean':clean}


setup(cmdclass=_commands,
      name='gemBS',
      version=__VERSION__,
      description='Python application to perform the different steps involved in the Bisulphite Pipeline.',
      author='Marcos Fernandez-Callejo, Santiago Marco-Sola, Simon Heath',
      author_email='marcos.fernandez@cnag.crg.eu',
      url='http://statgen.cnag.cat/gemBS/',
      packages=['gemBS'],
      package_data={"": [os.path.join("gemBS/gemBSbinaries", x) for x in ["readNameClean",
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

