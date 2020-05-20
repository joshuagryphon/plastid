#!/usr/bin/env python
"""Setup script for plastid. This is boilerplate except we needed to make
some workarounds to accomodate the fact that compiling C extensions requires
importing from cython, numpy, and pysam, which may not be pre-installed.

Specifically, the workaround:

- defines placeholder values for a number of variables (see [PLACEHOLDERS]),
  enabling `pip` to runthe `egg_info` command required for dependency parsing
  and installation. This allows `pip` to install cython, numpy, and pysam
  after a first pass through this file

- re-defines those placeholders if cython, numpy, and pysam are installed
  (see [REDEFINES]), allowing `pip` to build these extensions when it makes
  its second pass through the file

In addition:

- `build_ext`, `install`, and `develop` commands are wrapped to add a
  ``--recythonize`` flag, which causes Cython to regenerate ".c" files from
  ".pyx" files before the commands are run. This is merely a convenience for
  developers to force rebuild even if ".c" files appear up to date (but are not
  e.g. due to a switch in git repo)

- We added a `clean` command class was added to wipe old ".c" files,

"""
__author__ = "Joshua Griffin Dunn"

import os
import sys
import glob
from setuptools import setup, find_packages, Extension, Command
from setuptools.command.install import install
from setuptools.command.develop import develop
from pkg_resources import parse_version

plastid_version = "0.5.1"

# require python >= 2.7 (for 2.x) or >= 3.3 (for 3.x branch)
version_message = "plastid requires Python >= 2.7 or >= 3.3. Aborting installation."
ver = sys.version_info
if ver < (2, 7) or ver[0] == 3 and ver[1] < 3:
    raise RuntimeError(version_message)

#===============================================================================
# [CONSTANT]
#
# Setup that does not require pre-installation of C dependencies
#===============================================================================

# metadata ---------------------------------------------------------------------

with open("README.rst") as f:
    long_description = f.read()

setup_requires = [
    "numpy>=1.9.4",
    "pysam>=0.8.4",
    "cython>=0.22.0",
]

packages = find_packages() + [
    "kentUtils",
    "kentUtils.src.inc",
    "kentUtils.src.lib",
    "kentUtils.samtabix",
]

# trim dependencies if on readthedocs server, where many dependencies are mocked
on_rtd = os.environ.get("READTHEDOCS", None) == "True"

if not on_rtd:
    install_requires = [
        "scipy>=0.15.1",
        "pandas>=0.17.0",
        "matplotlib>=1.4.0",
        "biopython>=1.64",
        "twobitreader>=3.0.0",
        "termcolor",
    ] + setup_requires
else:
    install_requires = ["cython", "numpy", "pysam", "biopython", "termcolor"]


def get_scripts():
    """Detect command-line scripts automatically

    Returns
    -------
    list
        list of strings describing command-line scripts
    """
    binscripts = [X.replace(".py", "") for X in filter(lambda x: x.endswith(".py") and \
                                                                "__init__" not in x,
                                                      os.listdir(os.path.join("plastid",  "bin")))]
    return ["%s = plastid.bin.%s:main" % (X, X) for X in binscripts]


# required to build C extensions ----------------------------------------------

LIBRARIES = []
base_path = os.getcwd()

# embed method signatures and use Py2 or Py3 char specs, as appropriate
CYTHON_ARGS = {
    "embedsignature": True,
    "language_level": sys.version_info[0],
}

# we add a command class to force rebuild C files from pyx
# and a corresponding command-line argument
CYTHONIZE_COMMAND = "recythonize"  # command-class command from setup.py
CYTHONIZE_ARG = "recythonize"      # '--recythonize', passable to sdist, build_ext et c

# define options accepted by command-line argument
CYTHON_OPTIONS = [
        (CYTHONIZE_ARG,
         None,
         "If supplied,  use Cython to regenerate .c sources from pyx files " +\
         "(default False)"
        ),
]

# turn off --recythonize by default
CYTHON_DEFAULTS = [False]

# extension dependencies -------------------------------------------------------

# Build an extension of portions of Jim Kent's utilities,  which are needed by
# BigWigReader. The Kent utilities come under a permissive license that
# allows redistribution and modification for any use, including commercial,
# with the exception of src/portimpl.h, which appears to have its own license.
# Sources included here have therefore been modified not to depend on portimpl.h.
#
# JK's utilities are compiled into plastid.readers.bigwig, plastid.readers.bigbed,
# and plastid.readers.bbi_file

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

kent_deps  = [os.path.join(base_path, "kentUtils", "samtabix", X) for X in kent_samtabix]
kent_deps += [os.path.join(base_path, "kentUtils", "src", "lib", X) for X in kent_sources]

#===============================================================================
# [PLACEHOLDERS]
#
# Placeholders used to build egg_info if C dependencies are not installed.
#
# These placeholders enable `pip` to read this file, identify and install
# dependencies, and then re-run this file _after_ dependencies are installed to
# enable compilation of C extensions.
#
# If dependencies for compiled extensions are installed, these are overwritten
# or extended below in [REDEFINES]
#
#===============================================================================

# will be overwritten
ext_modules = []
DEFINE_MACROS = []
C_PATHS = []

# will be extended
INCLUDE_PATH = [
    os.path.join(base_path, "kentUtils", "src", "inc"),
    os.path.join(base_path, "kentUtils", "samtabix"),
]

command_classes = {}

#===============================================================================
# [REDEFINES]
#
# Detect if dependencies required to build C extensions are pre-installed
# if so, overwrite values from the section above
#===============================================================================

try:
    import numpy
    import pysam
    from Cython.Distutils import build_ext
    from Cython.Distutils.extension import Extension
    from Cython.Build import cythonize

    # extend include paths
    for mod in (numpy, pysam):
        ipart = mod.get_include()
        if isinstance(ipart, str):
            INCLUDE_PATH.append(ipart)
        elif isinstance(ipart, (list, tuple)):
            INCLUDE_PATH.extend(ipart)
        else:
            raise ValueError("Could not parse include path: %s" % ipart)

    DEFINE_MACROS = pysam.get_defines()

    # determine Pysam version
    # Several classes that we cimport moved in pysam 0.10.0,
    # so we define this macro in the build environment
    # to enable map_factories.pyx to cimport from the correct path
    CYTHON_COMPILE_TIME_ENV = {
        "PYSAM10": parse_version(pysam.__version__) >= parse_version("0.10.0")
    }

    # define extensions -------------------------------------------------------

    # redefine extensions
    extension_kwargs = {
        "include_dirs": INCLUDE_PATH,
        "language": "c",
        "cython_directives": CYTHON_ARGS,
    }

    # These extensions have no dependencies on kentUtils,
    # but do have dependencies on pysam
    noinclude_pyx = glob.glob(os.path.join(base_path, "plastid", "genomics", "*.pyx"))
    ext_modules = [
        Extension(
            x.replace(base_path + os.sep, "").replace(".pyx", "").replace(os.sep, "."), [x],
            libraries               = LIBRARIES,
            define_macros           = DEFINE_MACROS,
            cython_compile_time_env = CYTHON_COMPILE_TIME_ENV,
            **extension_kwargs
        ) for x in noinclude_pyx
    ] # yapf: disable

    # The following extensions do link to kentUtils, and also zlib
    bbifile = Extension(
        "plastid.readers.bbifile",
        ["plastid/readers/bbifile.pyx"] + kent_deps,
        libraries=LIBRARIES + ["z"],
        **extension_kwargs
    )

    bigwig = Extension(
        "plastid.readers.bigwig",
        ["plastid/readers/bigwig.pyx"] + kent_deps,
        libraries=LIBRARIES + ["z"],
        **extension_kwargs
    )

    bigbed = Extension(
        "plastid.readers.bigbed",
        ["plastid/readers/bigbed.pyx"] + kent_deps,
        libraries=LIBRARIES + ["z"],
        **extension_kwargs
    )

    ext_modules.append(bbifile)
    ext_modules.append(bigwig)
    ext_modules.append(bigbed)

    # define helper functions & classes for build -----------------------------

    # paths to sources
    PYX_PATHS = []
    for ex in ext_modules:
        PYX_PATHS.extend([X for X in ex.sources if X.endswith("pyx")])

    C_PATHS = [X.replace(".pyx", ".c") for X in PYX_PATHS]

    def wrap_command_classes(baseclass):
        """Add custom command-line `--recythonize` options to
        distutils/setuptools command classes

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
                for (op, _, _), default in zip(self.new_options, self.new_defaults):
                    setattr(self, op, default)

            def finalize_options(self):
                baseclass.finalize_options(self)
                if "--%s" % CYTHONIZE_ARG in self.distribution.script_args:
                    print("Cythonize on")
                    setattr(self, CYTHONIZE_ARG, True)

            def run(self):
                have_all = all([os.access(X, os.F_OK) for X in C_PATHS])
                if have_all == False:
                    print("Could not find .c files. Regenerating via recythonize.")
                    setattr(self, CYTHONIZE_ARG, True)
                    self.distribution.script_args.append("--%s" % CYTHONIZE_ARG)

                if getattr(self, CYTHONIZE_ARG) == True:
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
                if os.access(file_, os.F_OK):
                    print("clean_c_files: removing %s ..." % file_)
                    os.remove(file_)

    class build_c_from_pyx(build_ext):
        """Regenerate .c files from pyx files if --CYTHONIZE_ARG or
        CYTHONIZE_COMMAND is added to command line
        """
        user_options = build_ext.user_options + CYTHON_OPTIONS
        new_options = CYTHON_OPTIONS
        new_defaults = CYTHON_DEFAULTS

        cython_args = CYTHON_ARGS
        include_path = INCLUDE_PATH
        old_extensions = ext_modules

        description = "Regenerate .c files from .pyx source"

        def initialize_options(self):
            for (op, _, _), default in zip(self.new_options, self.new_defaults):
                setattr(self, op, default)
            build_ext.initialize_options(self)

        def finalize_options(self):
            if "--%s" % CYTHONIZE_ARG in self.distribution.script_args \
              or CYTHONIZE_COMMAND in self.distribution.script_args:

                self.run_command('clean')
                print("build_c_from_pyx: regenerating .c files from Cython")
                extensions = cythonize(
                    ext_modules,
                    compiler_directives=CYTHON_ARGS,
                    compile_time_env=CYTHON_COMPILE_TIME_ENV,
                )
                self.extensions = extensions

            build_ext.finalize_options(self)

        # override so that extensions are only built with build_ext, install, etc
        def run(self):
            pass

    # setup command classes ---------------------------------------------------

    command_classes = {
        CYTHONIZE_COMMAND: build_c_from_pyx,
        'build_ext': wrap_command_classes(build_ext),
        'install': wrap_command_classes(install),
        'develop': wrap_command_classes(develop),
        'clean': clean_c_files,
    }

except ImportError:
    print("plastid: Not all requirements pre-installed. Will need to bootstrap.")

#===============================================================================
# Program body
#===============================================================================

setup(

    name             = "plastid",
    version          = plastid_version,
    author           = "Joshua Griffin Dunn",
    author_email     = "joshua.g.dunn@gmail.com",
    maintainer       = "Joshua Griffin Dunn",
    maintainer_email = "joshua.g.dunn@gmail.com",
    long_description =  long_description,

    description      = "Tools for analysis of genomics & sequencing data",
    license          = "BSD 3-Clause",
    keywords         = "ribosome profiling riboseq rna-seq sequencing genomics biology",
    url              = "https://github.com/joshuagryphon/plastid",
    download_url     = "https://pypi.python.org/pypi/plastid/",
    platforms        = "OS Independent",

    classifiers      = [
         'Development Status :: 5 - Production/Stable',

         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3.3',
         'Programming Language :: Python :: 3.4',
         'Programming Language :: Python :: 3.5',
         'Programming Language :: Python :: 3.6',
         'Programming Language :: Python :: 3.7',
         'Programming Language :: Python :: 3.8',

         'Topic :: Scientific/Engineering :: Bio-Informatics',
         'Topic :: Software Development :: Libraries',

         'Intended Audience :: Science/Research',
         'Intended Audience :: Developers',

         'License :: OSI Approved :: BSD License',
         'Operating System :: POSIX',
         'Natural Language :: English',
    ],


    zip_safe = False,
    packages = packages,
    package_data = {
        "" : ["*.pyx", "*.pxd"],
        "kentUtils"          : ["README*"],
        "kentUtils.src.inc"  : ["*.h"],
        "kentUtils.src.lib"  : kent_sources + ["README"],
        "kentUtils.samtabix" : kent_samtabix + ["*.h", "COPYING", "README*", "AUTHORS"],
    },

    package_dir = {
        "plastid"            : "plastid",
        "kentUtils"          : "kentUtils",
        "kentUtils.src.inc"  : os.path.join("kentUtils", "src", "inc"),
        "kentUtils.src.lib"  : os.path.join("kentUtils", "src", "lib"),
        "kentUtils.samtabix" : os.path.join("kentUtils", "samtabix"),
    },

    entry_points = {
        "console_scripts" : get_scripts()
    },

    ext_modules      = ext_modules,
    cmdclass         = command_classes,

    setup_requires   = setup_requires,
    install_requires = install_requires,

    tests_require    = [ "nose>=1.0" ],
    test_suite       = "nose.collector",

) # yapf: disable
