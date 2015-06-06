#!/usr/bin/env python
from setuptools import setup, find_packages
import os
import yeti

with open("README.md") as f:
    long_description = f.read()


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



config_info = { "version"      : yeti.__version__,
                "entry_points" : dict(console_scripts=get_scripts()),
                "packages"     : find_packages(),
              }


setup(
    name = "yeti",
    install_requires = [
    	                "numpy>=1.9.2",
                        "scipy>=0.15.1",
                        "matplotlib>=1.3.0",
                        "biopython>=1.64",
                        "pysam>=0.7.7",
                        "pandas>=0.16.0",
                        ],
    include_package_data=True,
    package_data = {
#    "": ["*.*"],#
#     "": ["*.bam","*.bai","*.bed","*.bb","*.gtf","*.gff","*.gz","*.tbi","*.psl",
#          "*.bowtie","*.fa","*.wig","*.juncs","*.positions","*.sizes","*.as","*.txt"],
    },

    zip_safe = True,

    # metadata for upload to PyPI
    author           = "Joshua Griffin Dunn",
    author_email     = "joshua.g.dunn@gmail.com",
    maintainer       = "Joshua Griffin Dunn",
    maintainer_email = "Joshua Griffin Dunn",
    
    description = "Convert genomic datatypes into Pythonic objects useful to the SciPy stack",
    long_description = long_description,
    license   = "BSD 3-Clause",
    keywords  = "ribosome profiling riboseq rna-seq sequencing genomics biology",
    url       = "",   # project home page, if any
    platforms = "POSIX", # windows, 
    
    tests_require=["nose>=1.0"],
    test_suite = "nose.collector",
    
    classifiers=[
         'Development Status :: 4 - Beta',

         'Programming Language :: Python',
         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3.3',
         'Programming Language :: Python :: 3.4',

         'Topic :: Bioinformatics',
         'Topic :: Sequencing',
         'Topic :: Genomics'

         'Intended Audience :: End Users',
         'Intended Audience :: Biologists',
         'Intended Audience :: Computational Biologists',
         'Intended Audience :: Developers',

         'License :: BSD 3-Clause',
         
         'Operating System :: POSIX',
         #'Operating System :: MacOS :: MacOS X',
         #'Operating System :: Microsoft :: Windows',
        ],
    
    **config_info
)
