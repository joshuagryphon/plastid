#!/usr/bin/env python
"""Welcome to plastid!

This package contains various utilities for analyzing high-throughput sequencing
data, with an emphasis on simplicity for users. To this end, this package provides:

  #. A set of command-line scripts that implement common sequencing workflows
     (see |bin|).
  
  #. Readers that abstract data from various file formats into a minimal set of
     object types. These object types define APIs that easily interface with
     existing scientific tools, such as the `SciPy`_ stack (see |genomics| and
     |readers|)

  #. Tools to facilitate writing command-line scripts (see |scriptlib|)


Package overview
----------------
plastid is divided into the following subpackages:

    ==============    =========================================================
    Package           Contents
    --------------    ---------------------------------------------------------
    |bin|             Command-line scripts
    |genomics|        Classes and functions to manipulate genome annotations, alignments, and quantitative data
    |readers|         Parsers for various file formats
    |util|            Utilities (e.g. function decorators, exceptions, argument parsers)
    |test|            Unit and functional tests (requires download of test datasets)
    ==============    =========================================================
     
"""
__version__ = "0.3.2"
__author__  = "Joshua Griffin Dunn"
