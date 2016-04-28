#!/usr/bin/env python
"""Setup script for plastid. This is boilerplate except:

    1. Setup will error if numpy, Pysam, and Cython are not pre-installed.
       This is because this script needs those to run. This situation is
       is unavoidable at present given how package specification works.

    2. `build_ext`, `install`, and `develop` commands are wrapped to add
       a ``--recythonize`` flag, which causes Cython to regenerate ".c" files
       from ".pyx" files before the commands are run. This is merely a 
       convenience for developers to force rebuild even if ".c" files
       appear up to date (but are not e.g. due to a switch in git repo)

    3. a `clean` command class was added to wipe old ".c" files,
       e.g. optionally before meaking an sdist

    4. Kent utilities and dependencies are compiled and included for
       ``BigBedReader`` and ``BigWigReader``.

    5. Separately, if the package is being installed on readthedocs.org,
       it only installs absolutely necessary build dependencies, which are
       numpy and pysam.
"""
__author__ = "Joshua Griffin Dunn"

import os
import sys
import glob
import importlib
from setuptools import setup, find_packages, Extension, Command
from setuptools.command.install import install
from setuptools.command.develop import develop
from pkg_resources import parse_version

plastid_version = "0.4.6"

PYSAM_VER_NUM  = (0,8,4)
NUMPY_VER_NUM  = (1,9,0)
CYTHON_VER_NUM = (0,22,0)

PYSAM_VER_STR = parse_version(".".join([str(X) for X in PYSAM_VER_NUM]))
NUMPY_VER_STR = parse_version(".".join([str(X) for X in NUMPY_VER_NUM]))
CYTHON_VER_STR = parse_version(".".join([str(X) for X in CYTHON_VER_NUM]))

NUMPY_VERSION  = "numpy>=%s.%s.%s"  % NUMPY_VER_NUM
SCIPY_VERSION  = "scipy>=0.15.1"
CYTHON_VERSION = "cython>=%s.%s.%s" % CYTHON_VER_NUM
PYSAM_VERSION  = "pysam>=%s.%s.%s"  % PYSAM_VER_NUM

# require python >= 2.7 (for 2.x) or >= 3.3 (for 3.x branch)
version_message = "plastid requires Python >= 2.7 or >= 3.3. Aborting installation."
ver = sys.version_info
if ver < (2,7) or ver[0] == 3 and ver[1] < 3:
    raise RuntimeError(version_message)



# at present, setup/build can't proceed without numpy, Pysam, Cython
# adding these to setup_requires and/or install_requires doesn't work,
# because they are needed before then.
foundstr = "(found %s)"
nstr = pstr = cstr = "(no version found)"
numpyver = pysamver = cythonver = parse_version("0.0.0")
have_requirements = True

def check_version(modname):
    """Try to import a module by name.
    
    Parameters
    ----------
    modname : str
        Name of module

    Returns
    -------
    module or `None`
        imported module if it can be imported, otherwise `None`
    """
    try:
        mod = importlib.import_module(modname)
    except ImportError:
        mod = None
    return mod

try:
    numpy = check_version("numpy")
    if numpy is not None:
        numpyver = parse_version(numpy.__version__)
        nstr = foundstr % (numpy.__version__)

    pysam = check_version("pysam")
    if pysam is not None:
        pysamver = parse_version(pysam.__version__)
        pstr = foundstr % (pysam.__version__)

    Cython = check_version("Cython")
    if Cython is not None:
        cythonver = parse_version(Cython.__version__)
        cstr = foundstr % (Cython.__version__)
        
    if numpyver < NUMPY_VER_STR or pysamver < PYSAM_VER_STR or cythonver < CYTHON_VER_STR:
        raise ImportError()

    from Cython.Distutils import build_ext, Extension
    from Cython.Build import cythonize

except ImportError:
    ltmp = []
    ltmp.extend(NUMPY_VER_NUM)
    ltmp.append(nstr)
    ltmp.extend(PYSAM_VER_NUM)
    ltmp.append(pstr)
    ltmp.extend(CYTHON_VER_NUM)
    ltmp.append(cstr)

    print("""
    
*** IMPORTANT INSTALLATION INFORMATION ***

plastid setup requires the following to be preinstalled:

    numpy  %s.%s.%s or greater  %s
    Pysam  %s.%s.%s or greater  %s 
    Cython %s.%s.%s or greater  %s


Please install or upgrade these via pip, and retry:

    $ pip install --upgrade numpy pysam cython
    $ pip install plastid


If you believe you have more up-to-date versions than those detected above,
check the most up-to-date versions by typing:

    $ pip list | grep -E "(Cython|pysam|numpy)+"

If version numbers reported by `pip list` differ from those reported above,
please check your system to make sure there are not multiple versions 
of each package installed. This can be a problem e.g. if using anaconda.
In this case, Plastid can be installed and run within a virtual environment
(not a conda environment):

    # install virtualenv if you don't have it
    $ pip install virtualenv

    # make & activate virtualenv
    $ virtualenv /path/to/new/venv
    $ source /path/to/new/venv/bin/activate

    # install plastid in virtualenv
    $ pip install numpy pysam cython
    $ pip install --verbose plastid | tee install_log.txt


""" % tuple(ltmp))
    sys.exit(1)


#===============================================================================
# Simple stuff
#===============================================================================

with open("README.rst") as f:
    long_description = f.read()

setup_requires = [NUMPY_VERSION,PYSAM_VERSION,CYTHON_VERSION]
packages = find_packages()


def get_scripts():
    """Detect command-line scripts automatically

    Returns
    -------
    list
        list of strings describing command-line scripts
    """
    binscripts = [X.replace(".py","") for X in filter(lambda x: x.endswith(".py") and \
                                                                "__init__" not in x,
                                                      os.listdir(os.path.join("plastid","bin")))]
    return ["%s = plastid.bin.%s:main" % (X,X) for X in binscripts]



#===============================================================================
# Reconfigure build if on readthedocs.org
#
# If on doc server, only require essentials for building package:
#     Cython, numpy, pysam.
#
# Everything else can be mocked
#===============================================================================

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:
    install_requires = [
                        SCIPY_VERSION,
                        "pandas>=0.17.0",
                        "matplotlib>=1.4.0",
                        "biopython>=1.64",
                        "twobitreader>=3.0.0",
                        "termcolor",
                        ] + setup_requires
else:
    install_requires = ["cython","numpy","pysam","biopython","termcolor"]



#===============================================================================
# Paths to Cython extensions, included headers, et c
#===============================================================================

# options for compilation ------------------------------------------------------

CYTHONIZE_COMMAND = "recythonize" # command-class command from setup.py
CYTHONIZE_ARG = "recythonize"     # '--recythonize', passable to sdist, build_ext et c

CYTHON_ARGS= {
    "embedsignature" : True,
    "language_level" : sys.version_info[0],

    # enable/disable these for development/release
    #"linetrace" : True, # requires passing CYTHON_TRACE=1 to C compiler
    #"profile"   : True,
}


CYTHON_OPTIONS = [
        (CYTHONIZE_ARG,
         None,
         "If supplied, use Cython to regenerate .c sources from pyx files " +\
         "(default False)"
        ),
]

CYTHON_DEFAULTS = [False]
"""Default values for CYTHON_OPTIONS"""


# add headers from numpy & pysam to include paths for pyx, c files
base_path = os.getcwd()

INCLUDE_PATH = [os.path.join(base_path,"kentUtils","src","inc"),
                os.path.join(base_path,"kentUtils","samtabix")
                ]

for mod in (numpy,pysam):
    ipart = mod.get_include()
    if isinstance(ipart,str):
        INCLUDE_PATH.append(ipart)
    elif isinstance(ipart,(list,tuple)):
        INCLUDE_PATH.extend(ipart)
    else:
        raise ValueError("Could not parse include path: %s" % ipart)


# other flags; fill in as needed
LIBRARIES=[]
LIBRARY_DIRS=[]
RUNTIME_LIBRARY_DIRS=[]
EXTRA_OBJECTS=[]


#===============================================================================
# Build an extension of portions of Jim Kent's utilities, which are needed by
# BigWigReader. The Kent utilities come under a permissive license that
# allows redistribution and modification for any use, including commercial,
# with the exception of src/portimpl.h, which appears to have its own license.
# Sources included here have therefore been modified not to depend on portimpl.h.
#
# JK's utilities are compiled into plastid.readers.bigwig, plastid.readers.bigbed,
# and plastid.readers.bbi_file
#===============================================================================

kent_samtabix = [
    "kstring.c",
    "bgzf.c",
    "index.c",
]

kent_sources = [
    "aliType.c",
    "asParse.c",
    "base64.c",
    "basicBed.c",
    "bbiRead.c",
    "bbiWrite.c",
    "bigBed.c",
    "binRange.c",
    "bits.c",
    "bPlusTree.c",
    "bwgQuery.c",
    "bwgValsOnChrom.c",
    "cirTree.c",
    "common.c",
    "dlist.c",
    "dnautil.c",
    "dystring.c",
    "errAbort.c",
    "ffAli.c",
    "hash.c",
    "hmmstats.c",
    "https.c",
    "internet.c",
    "intExp.c",
    "kxTok.c",
    "linefile.c",
    "localmem.c",
    "memalloc.c",
    "obscure.c",
    "osunix.c",
    "pipeline.c",
    "psl.c",
    "rangeTree.c",
    "rbTree.c",
    "sqlNum.c",
    "sqlList.c",
    "tokenizer.c",
    "udc.c",
    "verbose.c",
    "wildcmp.c",
    "zlibFace.c",
]


kent_deps  = [os.path.join(base_path,"kentUtils","samtabix",X) for X in kent_samtabix]
kent_deps += [os.path.join(base_path,"kentUtils","src","lib",X) for X in kent_sources]


#===============================================================================
# Definition of extension modules
#===============================================================================


# Cython extensions without external dependencies
noinclude_pyx = glob.glob(os.path.join(base_path,"plastid","genomics","*.pyx"))
ext_modules = [Extension(x.replace(base_path+os.sep,"").replace(".pyx","").replace(os.sep,"."),
                         [x],
                         include_dirs=INCLUDE_PATH,
                         libraries=LIBRARIES,
                         library_dirs=LIBRARY_DIRS,
                         runtime_library_dirs=RUNTIME_LIBRARY_DIRS,
                         extra_objects=EXTRA_OBJECTS,
                         cython_directives=CYTHON_ARGS,
                        ) for x in noinclude_pyx]

bbifile = Extension(
    "plastid.readers.bbifile",
    ["plastid/readers/bbifile.pyx"] + kent_deps,
    language="c",
    include_dirs=INCLUDE_PATH,
    libraries=LIBRARIES + ["z"],
    library_dirs=LIBRARY_DIRS,
    runtime_library_dirs=RUNTIME_LIBRARY_DIRS,
    cython_directives=CYTHON_ARGS,
)

bigwig = Extension(
    "plastid.readers.bigwig",
    ["plastid/readers/bigwig.pyx"] + kent_deps,
    language="c",
    include_dirs=INCLUDE_PATH,
    libraries=LIBRARIES + ["z"],
    library_dirs=LIBRARY_DIRS,
    runtime_library_dirs=RUNTIME_LIBRARY_DIRS,
    cython_directives=CYTHON_ARGS,
)

bigbed = Extension(
    "plastid.readers.bigbed",
    ["plastid/readers/bigbed.pyx"] + kent_deps,
    language="c",
    include_dirs=INCLUDE_PATH,
    libraries=LIBRARIES + ["z"],
    library_dirs=LIBRARY_DIRS,
    runtime_library_dirs=RUNTIME_LIBRARY_DIRS,
    cython_directives=CYTHON_ARGS,
)

ext_modules.append(bbifile)
ext_modules.append(bigwig)
ext_modules.append(bigbed)


# paths to sources
PYX_PATHS = []
for ex in ext_modules:
    PYX_PATHS.extend([X for X in ex.sources if X.endswith("pyx")])

C_PATHS   = [X.replace(".pyx",".c") for X in PYX_PATHS]


# classes & functions for compilation -----------------------------------------

def wrap_command_classes(baseclass):
    """Add custom command-line `--recythonize` options to distutils/setuptools
    command classes
    
    Parameters
    ----------
    baseclass : setuptools.Command or subclass (not instance!)
        Command class to modify

    Returns
    -------
    class
        Modified class
    """
    class subclass(baseclass):
        user_options = baseclass.user_options + CYTHON_OPTIONS
        new_options = CYTHON_OPTIONS
        new_defaults = CYTHON_DEFAULTS

        def initialize_options(self):
            baseclass.initialize_options(self)
            for (op,_,_), default in zip(self.new_options,
                    self.new_defaults):
                setattr(self,op,default)
        
        def finalize_options(self):
            baseclass.finalize_options(self)
            if "--%s" % CYTHONIZE_ARG in self.distribution.script_args:
                print("Cythonize on")
                setattr(self,CYTHONIZE_ARG,True)

        def run(self):
            have_all = all([os.access(X,os.F_OK) for X in C_PATHS])
            if have_all == False:
                print("Could not find .c files. Regenerating via recythonize.")
                setattr(self,CYTHONIZE_ARG,True)
                self.distribution.script_args.append("--%s" % CYTHONIZE_ARG)
            
            if getattr(self,CYTHONIZE_ARG) == True:
                print("Cythonizing")
                self.run_command(CYTHONIZE_COMMAND)

            return baseclass.run(self)

    subclass.__name__ = "cython_%s" % baseclass.__name__
    return subclass


class clean_c_files(Command):
    """Remove previously generated .c files"""
    user_options = []
    description = "Remove previously generated .c files"
 
    def initialize_options(self):
        pass
 
    def finalize_options(self):
        pass
 
    def run(self):
        for file_ in C_PATHS:
            if os.access(file_,os.F_OK):
                print("clean_c_files: removing %s ..." % file_)
                os.remove(file_)

 
class build_c_from_pyx(build_ext):
    """Regenerate .c files from pyx files if --CYTHONIZE_ARG or CYTHONIZE_COMMAND is added to command line"""
    user_options = build_ext.user_options + CYTHON_OPTIONS
    new_options  = CYTHON_OPTIONS
    new_defaults = CYTHON_DEFAULTS

    cython_args  = CYTHON_ARGS
    include_path = INCLUDE_PATH
    old_extensions = ext_modules

    description = "Regenerate .c files from .pyx source"

    def initialize_options(self):
        for (op,_,_), default in zip(self.new_options,
                self.new_defaults):
            setattr(self,op,default)
        build_ext.initialize_options(self)

    def finalize_options(self):
        if "--%s" % CYTHONIZE_ARG in self.distribution.script_args or CYTHONIZE_COMMAND in self.distribution.script_args:
            self.run_command('clean')
            print("build_c_from_pyx: regenerating .c files from Cython")
            extensions = cythonize(ext_modules,compiler_directives=CYTHON_ARGS)
            self.extensions = extensions

        build_ext.finalize_options(self)

    def run(self): # override so that extensions are only built with build_ext, install, etc
        pass



#===============================================================================
# Program body 
#===============================================================================

FAIL_MSG = \
"""Setup failed. If this is due to a compiler error, try rebuilding the Cython
sources. If installing via pip, make sure all dependencies are already
pre-installed, then separately run:

    $ pip install --verbose --upgrade --install-option="--recythonize" plastid | tee install_log.txt


Or, if installing manually via setup.py:

    $ python setup.py install --recythonize | tee install_log.txt
"""

setup(

    # package metadata
    name             = "plastid",
    version          = plastid_version,
    author           = "Joshua Griffin Dunn",
    author_email     = "joshua.g.dunn@gmail.com",
    maintainer       = "Joshua Griffin Dunn",
    maintainer_email = "Joshua Griffin Dunn",
    long_description =  long_description,

    description      = "Tools for analysis of genomics & sequencing data",
    license          = "BSD 3-Clause",
    keywords         = "ribosome profiling riboseq rna-seq sequencing genomics biology",
    url              = "https://github.com/joshuagryphon/plastid",
    download_url     = "https://pypi.python.org/pypi/plastid/",
    platforms        = "OS Independent",
 
    classifiers      = [
         'Development Status :: 4 - Beta',

         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3.3',
         'Programming Language :: Python :: 3.4',
         'Programming Language :: Python :: 3.5',

         'Topic :: Scientific/Engineering :: Bio-Informatics',
         'Topic :: Software Development :: Libraries',

         'Intended Audience :: Science/Research',
         'Intended Audience :: Developers',

         'License :: OSI Approved :: BSD License',
         'Operating System :: POSIX',
         'Natural Language :: English',
    ],


    # packaging info
    zip_safe = False,
    packages = find_packages() + ["kentUtils",
                                  "kentUtils.src.inc",
                                  "kentUtils.src.lib",
                                  "kentUtils.samtabix"],

    package_data = { 
        "" : ["*.pyx","*.pxd"],
        "kentUtils"          : ["README*"],
        "kentUtils.src.inc"  : ["*.h"],
        "kentUtils.src.lib"  : kent_sources + ["README"],
        "kentUtils.samtabix" : kent_samtabix + ["*.h","COPYING","README*","AUTHORS"],
    },

    package_dir = {
        "plastid"            : "plastid",
        "kentUtils"          : "kentUtils",
        "kentUtils.src.inc"  : os.path.join("kentUtils","src","inc"),
        "kentUtils.src.lib"  : os.path.join("kentUtils","src","lib"),
        "kentUtils.samtabix" : os.path.join("kentUtils","samtabix"),
    },
 
    entry_points = {
        "console_scripts" : get_scripts()
    },

    ext_modules = ext_modules,
    
    cmdclass    = {
        CYTHONIZE_COMMAND : build_c_from_pyx,
        'build_ext' : wrap_command_classes(build_ext), 
        'install'   : wrap_command_classes(install),
        'develop'   : wrap_command_classes(develop),
        'clean'     : clean_c_files,
    },

    setup_requires   = setup_requires,
    install_requires = install_requires,
   
    tests_require     = [ "nose>=1.0" ],
    test_suite        = "nose.collector",

)
