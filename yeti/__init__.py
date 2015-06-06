#!/usr/bin/env python
"""Welcome to yeti!

This package contains various utilities for analyzing high-throughput sequencing
data, with an emphasis on simplicity for users. To this end, this package provides:

  #. A set of command-line scripts that implement common sequencing workflows
     (see |bin|).
  
  #. Readers that abstract data from various file formats into a minimal set of
     object types. These object types define APIs that easily interface with
     existing scientific tools, such as the `SciPy`_ stack
     (see |readers|), facilitating analysis.

  #. Tools for quickly writing command-line scripts that leverage the object
     types defined in this package.


Package overview
----------------
yeti is divided into the following sub-packages:

    ==============    =========================================================
    Package           Contents
    --------------    ---------------------------------------------------------
    |bin|             Command-line scripts
    |genomics|        Object types that model genome annotations and quantitative data
    |readers|         Parsers for various file formats that yield objects of the types found in |genomics|
    |util|            Utilities (e.g. function decorators, exceptions, argument parsers)
    |test|            Unit and functional tests
    ==============    =========================================================
     
"""
__version__ = "0.1.0"
__author__  = "Joshua Griffin Dunn"
