#!/usr/bin/env python
import os
import sys
import glob
from collections import Iterable
from setuptools import setup, find_packages, Extension, Command
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install

NUMPY_VERSION  = "numpy>=1.9"
SCIPY_VERSION  = "scipy>=0.15.1"
CYTHON_VERSION = "cython>=0.22"
PYSAM_VERSION  = "pysam>=0.8.3"

# at present, setup/build can't proceed without Cython, numpy, & Pysam
# adding these to setup_requires and/or install_requires doesn't work,
# because they are needed before them.
try:
#    from Cython.Build import cythonize
#    from Cython.Distutils import build_ext
    import numpy
    import pysam
except ImportError:
    print("""plastid setup requires %s, %s, and %s to be installed. Please
install these via pip, and retry:

    $ pip install cython numpy pysam
""" % (CYTHON_VERSION,NUMPY_VERSION,PYSAM_VERSION))
    sys.exit(1)


#===============================================================================
# Constants/variables used below
#===============================================================================

use_cython = False

setup_requires = [NUMPY_VERSION,PYSAM_VERSION]
packages = find_packages()

with open("README.rst") as f:
    long_description = f.read()


base_path = os.getcwd()


pyx_paths = [os.path.join(base_path,"plastid","genomics","*.pyx")]
"""Paths to Cython .pyx files"""

c_paths = glob.glob(os.path.join(base_path,"plastid","genomics","*.c"))
"""Potential path to cythonized .c files"""

# append numpy & pysam headers to include path for gcc & Cython compilers
ipath = []
for mod in (numpy,pysam):
    ipart = mod.get_include()
    if isinstance(ipart,str):
        ipath.append(ipart)
    elif isinstance(ipart,(list,tuple)):
        ipath.extend(ipart)
    else:
        raise ValueError("Could not parse include path: %s" % ipart)

# define extensions for .c files
ext_modules = [Extension(X.replace(base_path+os.sep,"").replace(".c","").replace(os.sep,"."),
                        [X],include_dirs=ipath) for X in c_paths]

#ext_modules = [Extension("plastid.genomics.c_common",[base+"/plastid/genomics/c_common.c"],include_dirs=ipath),
#Extension("plastid.genomics.roitools",[base+"/plastid/genomics/roitools.c"],include_dirs=ipath),
#Extension("plastid.genomics.map_factories",[base+"/plastid/genomics/map_factories.c"],include_dirs=ipath),
#]

# if using Cython to regenerate c files:
if use_cython == True:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
    import numpy
    import pysam

    # required for setup.py to run
    setup_requires = [CYTHON_VERSION,
                      PYSAM_VERSION,
                      NUMPY_VERSION]

    cython_args= {
        "embedsignature" : True,
        "language_level" : sys.version_info[0],

        # enable/disable these for development/release
        #"linetrace" : True, # requires passing CYTHON_TRACE=1 to C compiler
        #"profile"   : True,
    }
    """Cython compiler options"""

    ext_modules = cythonize([Extension("*",pyx_paths,include_dirs = ipath)],
                             compiler_directives = cython_args,
                             include_path = ipath)        

# if on doc server, only require Cython, numpy, pysam. Everything else
# is simulated via mock
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:
     install_requires = [
                        SCIPY_VERSION,
                        "pandas>=0.16.0",
                        "matplotlib>=1.3.0",
                        "biopython>=1.64",
                        "twobitreader>=3.0.0",
                        ] + setup_requires
else:
    install_requires = ["cython","numpy","pysam"]


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

class install_after_build(install):
    def run(self):
        self.run_command("build_ext")

        return setuptools.command.install.install.run(self)

setup(

    # package metadata
    name             = "plastid",
    version          = "0.3.0", #plastid.__version__,
    author           = "Joshua Griffin Dunn",
    author_email     = "joshua.g.dunn@gmail.com",
    maintainer       = "Joshua Griffin Dunn",
    maintainer_email = "Joshua Griffin Dunn",
    long_description =  long_description,

    description      = "Convert genomic datatypes into Pythonic objects useful to the SciPy stack",
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
    #include_package_data = True,
    package_data         = { "" : ["*.pyx","*.pxd","*.c"], },
#    exclude_package_data = { "" : ["*.c"], },
#    "": ["*.*"],#
#     "": ["*.bam","*.bai","*.bed","*.bb","*.gtf","*.gff","*.gz","*.tbi","*.psl",
#          "*.bowtie","*.fa","*.wig","*.juncs","*.positions","*.sizes","*.as","*.txt"],

    # dependencies, entrypoints, ext modules, et c
    entry_points     = { "console_scripts" : get_scripts()
                       },

    ext_modules = ext_modules,
    cmdclass    = {
                    #'build_ext' : CythonBuildExtProxy,
                    'build_ext' : build_ext,
                    #'install' : install_after_build,
                   },

    setup_requires   = setup_requires,
    install_requires = install_requires,
   
    tests_require     = [ "nose>=1.0" ],
    test_suite        = "nose.collector",

)
