#from setuptools import setup, find_packages
import os
import yeti
import sys
import numpy
import pysam

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

#===============================================================================
# Constants/variables used below
#===============================================================================

NUMPY_VERSION  = "numpy>=1.9.2"
SCIPY_VERSION  = "scipy>=0.15.1"
CYTHON_VERSION ="cython>=0.22"

include_path = []
"""Paths to headers/sources used by pyx files"""

for mod in (numpy,pysam):
    ic = mod.get_include()
    if isinstance(ic,list):
        include_path.extend(ic)
    elif isinstance(ic,str):
        include_path.append(ic)
    else:
        raise TypeError("Includes for %s are '%s' (type %s)" % (mod.__name__,ic,type(ic)))


ext_modules = [
    Extension("*",["yeti/genomics/*.pyx"],include_dirs=include_path)
]
"""Cython modules to compile"""


cython_directives = {
    "embedsignature" : True,
    "language_level" : sys.version_info[0],

    # enable/disable these for development/release
    "linetrace" : True, # requires passing CYTHON_TRACE=1 to C compiler
    "profile"   : True,
}
"""Cython compiler options"""


with open("README.rst") as f:
    long_description = f.read()


#===============================================================================
# Helper functions for setup
#===============================================================================

def get_scripts():
    """Detect command-line scripts before egg build
    
    Returns
    -------
    list
        list of strings describing command-line scripts
    """
    binscripts = [X.replace(".py","") for X in filter(lambda x: x.endswith(".py"),
                                                      os.listdir(os.path.join("yeti","bin")))]
    return ["%s = yeti.bin.%s:main" % (X,X) for X in binscripts]


#===============================================================================
# Program body 
#===============================================================================

setup(

    # package metadata
    name             = "yeti",
    version          = yeti.__version__,
    author           = "Joshua Griffin Dunn",
    author_email     = "joshua.g.dunn@gmail.com",
    maintainer       = "Joshua Griffin Dunn",
    maintainer_email = "Joshua Griffin Dunn",
    long_description =  long_description,

    description      = "Convert genomic datatypes into Pythonic objects useful to the SciPy stack",
    license          = "BSD 3-Clause",
    keywords         = "ribosome profiling riboseq rna-seq sequencing genomics biology",
    url              = "", # github link
    download_url     = "", # PyPI link
    platforms        = "OS Independent",
 
    classifiers      = [
         'Development Status :: 4 - Beta',

         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3.3',
         'Programming Language :: Python :: 3.4',

         'Topic :: Scientific/Engineering :: Bio-Informatics',
         'Topic :: Software Development :: Libraries',

         'Intended Audience :: Intended Audience :: Science/Research',
         'Intended Audience :: Developers',

         'License :: OSI Approved :: BSD License',
         #'Operating System :: OS Independent',
         'Operating System :: POSIX',
         'Natural Language :: English',
    ],


    # packaging info
    zip_safe             = True,
    packages             = ["yeti"], # find_packages()
    include_package_data = True,
    package_data         = {
#    "": ["*.*"],#
#     "": ["*.bam","*.bai","*.bed","*.bb","*.gtf","*.gff","*.gz","*.tbi","*.psl",
#          "*.bowtie","*.fa","*.wig","*.juncs","*.positions","*.sizes","*.as","*.txt"],
                           },


    # dependencies, entrypoints, ext modules, et c
    cmd_class        = {
                         'build_ext' : build_ext,
                       },

    entry_points     = { "console_scripts" : get_scripts()
                       },

    ext_modules      = cythonize(ext_modules,
                                 include_path=include_path,
                                 compiler_directives=cython_directives,
                                ),

    install_requires = [
    	                NUMPY_VERSION, #"numpy>=1.9.2",
                        SCIPY_VERSION, #"scipy>=0.15.1",
                        "matplotlib>=1.3.0",
                        "biopython>=1.64",
                        "pysam>=0.7.7",
                        "pandas>=0.16.0",
                        "twobitreader>=3.0.0",
                        CYTHON_VERSION, #"cython>=0.22",
                        ],
   
    tests_require     = [ "nose>=1.0" ],
    test_suite        = "nose.collector",

)
