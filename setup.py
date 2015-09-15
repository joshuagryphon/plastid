#!/usr/bin/env python
import os
import plastid
import sys
import glob
from collections import Iterable

from distutils.command.build_ext import build_ext as d_build_ext
from setuptools import setup, find_packages, Extension, Command

#===============================================================================
# Constants/variables used below
#===============================================================================


on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
packages = find_packages()

with open("README.rst") as f:
    long_description = f.read()


NUMPY_VERSION  = "numpy>=1.9"
SCIPY_VERSION  = "scipy>=0.15.1"
CYTHON_VERSION = "cython>=0.22"
PYSAM_VERSION  = "pysam>=0.8.3"

cython_args= {
    "embedsignature" : True,
    "language_level" : sys.version_info[0],

    # enable/disable these for development/release
    "linetrace" : True, # requires passing CYTHON_TRACE=1 to C compiler
    "profile"   : True,
}
"""Cython compiler options"""

ext_paths = ["plastid/genomics/*.pyx"]
"""Paths to Cython .pyx files"""

c_paths = [os.path.join(X.replace(",",os.sep),"*.c") for X in packages]
"""Potential paths to cythonized .c files files"""

ipath = []
setup_requires = []
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy
import pysam
setup_requires = [CYTHON_VERSION,
                  PYSAM_VERSION,
                  NUMPY_VERSION]
ipath = []
for mod in (numpy,pysam):
    ipart = mod.get_include()
    if isinstance(ipart,str):
        ipath.append(ipart)
    elif isinstance(ipart,Iterable):
        ipath.extend(ipart)
    else:
        raise ValueError("Could not parse include path: %s" % ipart)

    ext_modules = cythonize([Extension("*",ext_paths,include_dirs = ipath)],
                             compiler_directives = cython_args,
                             include_path = ipath)        
if not on_rtd:
     install_requires = [
                        SCIPY_VERSION,
                        "pandas>=0.16.0",
                        "matplotlib>=1.3.0",
                        "biopython>=1.64",
                        "twobitreader>=3.0.0",
                        ] + setup_requires
#    ipath = []
#    for mod in (numpy,pysam):
#        ipart = mod.get_include()
#        if isinstance(ipart,str):
#            ipath.append(ipart)
#        elif isinstance(ipart,Iterable):
#            ipath.extend(ipart)
#        else:
#            raise ValueError("Could not parse include path: %s" % ipart)
#
#        ext_modules = cythonize([Extension("*",ext_paths,include_dirs = ipath)],
#                                 compiler_directives = cython_args,
#                                 include_path = ipath)        
else:
#    from setuptools.command import build_ext
#    ext_modules = []
    install_requires = ["cython","numpy","pysam"]
                      #  "scipy","numpy","pandas","matplotlib","biopython","#  ]


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
    binscripts = [X.replace(".py","") for X in filter(lambda x: x.endswith(".py") and \
                                                                "__init__" not in x,
                                                      os.listdir(os.path.join("plastid","bin")))]
    return ["%s = plastid.bin.%s:main" % (X,X) for X in binscripts]
   

#class CythonBuildExtProxy(d_build_ext):
#    """Proxy class to enable use of Cython's `build_ext` command class without
#    having to import it before Cython is installed
#    """
#    def __init__(self,*args,**kwargs):
#        # import and pull a switcheroo on type
#        from Cython.Distutils import build_ext
#        self._obj = build_ext(*args,**kwargs)
#
#    def __str__(self):
#        return str(self._obj)
#
#    def __getattr__(self,attr):
#        return getattr(self._obj,attr)
#
#    def __setattr__(self,attr,val):
#        return setattr(self._obj,attr,val)
#
#class LateCython(list):
#    """Workaround:
#    
#    numpy, pysam, and Cython all need to be installed for the .pyx files
#    to build. Putting these in the `setup()` argument `setup_requires`
#    is insufficient, because other arguments passed to `setup()`
#    require numpy, pysam, and Cython, which won't yet be installed.
#
#    So, this class delays evaluation of those arguments until all
#    dependencies in `setup_requires` have been installed, allowing
#    allowing numpy, pysam, and Cython to be imported, and the pyx files
#    to be built.
#
#    Idea for proxy class from:
#        https://stackoverflow.com/questions/11010151/distributing-a-shared-library-and-some-c-code-with-a-cython-extension-module/26698408#26698408
#    """
#    def __init__(self,ext_paths,cython_kwargs={}):
#        """
#        Parameters
#        ----------
#        ext_paths : list
#            List of paths to Cython .pyx files
#
#        cython_kwargs : dict, optional
#            Cython compiler directives
#        """
#        self.ext_paths = ext_paths
#        self.cython_kwargs = cython_kwargs
#        self._extensions = None
#
#    def __getattr__(self,attr):
#        if attr == "extensions":
#            if self._extensions is None:
#                self._cythonize()
#            return self._extensions
#
#    def _cythonize(self):
#        """Add packages from `setup_requires` to sys.path, and then
#        cythonize extensions. Only called when contents of `self` are iterated,
#        i.e. when setuptools is building the package extensions.
#        """
#        if self._extensions is None:
#            try:
#            # add all eggs from `setup_requires` to sys.path
#                eggs = glob.glob(".eggs/*egg")
#                for n,egg in enumerate(eggs):
#                    print("Adding %s to sys.path" % egg)
#                    sys.path.append(egg)
#
#                # get include paths
#                import numpy
#                import pysam
#                ipath = []
#                for mod in (numpy,pysam):
#                    ipart = mod.get_include()
#                    if isinstance(ipart,str):
#                        ipath.append(ipart)
#                    elif isinstance(ipart,Iterable):
#                        ipath.extend(ipart)
#                    else:
#                        raise ValueError("Could not parse include path: %s" % ipart)
#
#                # try to build if not yet built
#                from Cython.Build import cythonize
#                extensions = [Extension("*",self.ext_paths,include_dirs = ipath)]
#                b = cythonize(extensions,
#                              include_path = ipath,
#                              compiler_directives = self.cython_kwargs)
#                self._extensions = b
#
#                for i in range(n+1):
#                    sys.path.pop()
#                return self._extensions
#            except ImportError:
#                for i in range(n+1):
#                    sys.path.pop()
#
#        return self._extensions
#
#    def __iter__(self):
#        return iter(self.extensions)
#        
#    def __len__(self):
#        return len(self.ext_paths)
#
#    def __getitem__(self,key):
#        return self.extensions[k]


#===============================================================================
# Program body 
#===============================================================================


setup(

    # package metadata
    name             = "plastid",
    version          = plastid.__version__,
    author           = "Joshua Griffin Dunn",
    author_email     = "joshua.g.dunn@gmail.com",
    maintainer       = "Joshua Griffin Dunn",
    maintainer_email = "Joshua Griffin Dunn",
    long_description =  long_description,

    description      = "Convert genomic datatypes into Pythonic objects useful to the SciPy stack",
    license          = "BSD 3-Clause",
    keywords         = "ribosome profiling riboseq rna-seq sequencing genomics biology",
    url              = "https://github.com/joshuagryphon/plastid", # github link
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
    packages             = packages,
    include_package_data = True,
    package_data         = {
#    "": ["*.*"],#
#     "": ["*.bam","*.bai","*.bed","*.bb","*.gtf","*.gff","*.gz","*.tbi","*.psl",
#          "*.bowtie","*.fa","*.wig","*.juncs","*.positions","*.sizes","*.as","*.txt"],
                           },

    # dependencies, entrypoints, ext modules, et c
    entry_points     = { "console_scripts" : get_scripts()
                       },

    # the following defer import of Cython, numpy, and pysam
    # until after setup() has processed `setup_requires`

    #ext_modules = LateCython(ext_paths,cython_args),
    ext_modules = ext_modules,
    cmdclass    = {
                    #'build_ext' : CythonBuildExtProxy,
                    'build_ext' : build_ext,
                   },


    setup_requires   = setup_requires,
    install_requires = install_requires,
   
    tests_require     = [ "nose>=1.0" ],
    test_suite        = "nose.collector",

)
