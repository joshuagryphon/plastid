#!/usr/bin/env python
import os
import yeti
import sys
from collections import Iterable

from distutils.command.build_ext import build_ext as d_build_ext
from setuptools import setup, find_packages, Extension, Command


#===============================================================================
# Constants/variables used below
#===============================================================================

NUMPY_VERSION  = "numpy>=1.9.2"
SCIPY_VERSION  = "scipy>=0.15.1"
CYTHON_VERSION ="cython>=0.22"

with open("README.rst") as f:
    long_description = f.read()

ext_paths = ["yeti/genomics/*.pyx"]
"""Paths to Cython .pyx files"""

cython_args= {
    "embedsignature" : True,
    "language_level" : sys.version_info[0],

    # enable/disable these for development/release
    "linetrace" : True, # requires passing CYTHON_TRACE=1 to C compiler
    "profile"   : True,
}
"""Cython compiler options"""


#===============================================================================
# Helper functions/classes for setup
#===============================================================================

def get_scripts():
    """Detect command-line scripts automatically

    Returns
    -------
    list
        list of strings describing command-line scripts
    """
    binscripts = [X.replace(".py","") for X in filter(lambda x: x.endswith(".py"),
                                                      os.listdir(os.path.join("yeti","bin")))]
    return ["%s = yeti.bin.%s:main" % (X,X) for X in binscripts]
    

class CythonBuildExtProxy(Command):
    """Proxy class to enable use of Cython's `build_ext` without
    having to import it before Cython is installed
    """
    user_options = d_build_ext.user_options

    def __new__(cls):
        obj = super(BuildExtClass,cls).__new__(cls)
        return obj

    def __init__(self,*args,**kwargs):
        from Cython.Distutils import build_ext
        self.__class__ = build_ext
        build_ext.__init__(self,*args,**kwargs)


class LateCython(list):
    """Workaround:
    
    numpy, pysam, and Cython all need to be installed for the .pyx files
    to build. Putting these in the setup() argument `setup_requires`
    is insufficient, because other arguments passed to `setup_requires`
    require numpy, pysam, and Cython.

    So, this class delays evaluation of those arguments until all
    dependencies in `setup_requires` have been installed, allowing
    allowing numpy, pysam, and Cython to be imported.

    Idea for proxy class from:
        https://stackoverflow.com/questions/11010151/distributing-a-shared-library-and-some-c-code-with-a-cython-extension-module/26698408#26698408
    """
    def __init__(self,ext_paths,cython_kwargs={}):
        """
        Parameters
        ----------
        ext_paths : list
            List of paths to Cython .pyx files

        cython_kwargs : dict, optional
            Cython compiler directives
        """
        from Cython.Distutils import build_ext
        self.ext_paths = ext_paths
        self.cython_kwargs = cython_kwargs
        self.build_ext = build_ext

    def get_include_path(self):
        """Get path to headers & pyx files for dependencies

        Returns
        -------
        list
            List of paths to pass to cython compiler
        """
        import numpy
        import pysam
        ipath = []
        for mod in (numpy,pysam):
            ipart = mod.get_include()
            if isinstance(ipart,str):
                ipath.append(ipart)
            elif isinstance(ipart,Iterable):
                ipath.extend(ipart)
            else:
                raise ValueError("Could not parse include path: %s" % ipart)

        return ipath

    def _construct(self):
        """Construct list of extensions and dependencies. Only called
        when contents of `self` are iterated.
        """
        try:
            from Cython.Build import cythonize
            include_path = self.get_include_path()
            extensions = [Extension("*",self.ext_paths,include_dirs=include_path)]
            b = cythonize(extensions,
                          include_path = include_path,
                          compiler_directives = self.cython_kwargs)

            return b
        except ImportError:
            raise RuntimeError("Cython not installed. Please install %s" % CYTHON_VERSION)

    def __iter__(self):
        return iter(self._construct())
        
    def __len__(self):
        return 1

    def __getitem__(self,key):
        return self._construct()


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
    packages             = find_packages()
    include_package_data = True,
    package_data         = {
#    "": ["*.*"],#
#     "": ["*.bam","*.bai","*.bed","*.bb","*.gtf","*.gff","*.gz","*.tbi","*.psl",
#          "*.bowtie","*.fa","*.wig","*.juncs","*.positions","*.sizes","*.as","*.txt"],
                           },

    # dependencies, entrypoints, ext modules, et c
    entry_points     = { "console_scripts" : get_scripts()
                       },

    # these classes defer import of Cython, numpy, and pysam
    # until after setup() has installed them, processing
    # `setup_requires` and `install_requires`
    ext_modules = LateCython(ext_paths,cython_args),
    cmdclass    = {
                    'build_ext' : CythonBuildExtProxy,
                   },


    setup_requires   = [CYTHON_VERSION,
                        NUMPY_VERSION,
                        SCIPY_VERSION],

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
