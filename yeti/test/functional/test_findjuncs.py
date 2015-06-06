#!/usr/bin/env python
"""Test suite for :py:mod:`yeti.bin.findjuncs`"""
import tempfile
import os

from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from yeti.test.functional.base import execute_helper
from yeti.test.ref_files import REF_FILES
from yeti.bin.findjuncs import main
from yeti.util.services.decorators import catch_stderr
from yeti.util.services.mini2to3 import cStringIO


#===============================================================================
# INDEX: global constants used by tests
#===============================================================================

test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "yeti.bin.findjuncs",
    "ref_file_path"  : resource_filename("yeti","test/data/"),
    "temp_file_path" : tempfile.mkdtemp(prefix="findjuncs"),
}

# Define tests as tuples of:
# 1. Command-line style arguments to pass to ``main``
# 2. A list of reference files that output should be compared against
# 3. A list of output files created by running ``main`` with the arguments provided in (1)
# 4. A list of strings specifying how equality should be evaluated
basename = os.path.join(test_info["temp_file_path"],"test_findjuncs")

findjuncs_tests = [
    # find all junctions in one GTF
    ("%s --annotation_files %s --annotation_format GTF2" % (basename,REF_FILES["yeast_gtf"]),
     [REF_FILES["yeast_findjuncs_bed"]],
     ["%s.bed" % basename],
     ["--no_header --sort_keys 3"]),

    # find all junctions in two GTFs; report only once
    ("%s --annotation_files %s %s" % (basename,REF_FILES["yeast_gtf"],REF_FILES["yeast_gtf"]),
     [REF_FILES["yeast_findjuncs_bed"]],
     ["%s.bed" % basename],
     ["--no_header --sort_keys 3"]),

    # find all junctions in one GTF, export tophat juncs
    ("%s --annotation_files %s --export_tophat" % (basename,REF_FILES["yeast_gtf"]),
     [REF_FILES["yeast_findjuncs_bed"],
      REF_FILES["yeast_findjuncs_juncs"],
      ],
     ["%s.bed" % basename,
      "%s.juncs" % basename],
     ["--no_header --sort_keys 3",
      "--no_header"]),
    
    # pool different but overlapping sets of junctions from two BED files, export tophat juncs
    ("%s --annotation_files %s %s --annotation_format BED --export_tophat" % (basename,
                                                                              REF_FILES["yeast_findjuncs_top"],
                                                                              REF_FILES["yeast_findjuncs_bot"]),
      [REF_FILES["yeast_findjuncs_bed"],
       REF_FILES["yeast_findjuncs_juncs"]],
     ["%s.bed" % basename,
      "%s.juncs" % basename],
     ["--no_header --sort_keys 3",
      "--no_header"]),
]
"""Functional tests of :py:mod:`yeti.bin.findjuncs`.

Tests are specified as tuples of:

    1. Command-line style arguments to pass to :py:func:`main`

    2. A list of reference files that output should be compared against

    3. A list of output files created by running :py:func:`main`
       with the arguments provided in (1)

    4. A list of strings specifying how equality should be evaluated
"""

@attr(test="functional")
@attr(speed="slow")
def do_test():
    """Perform functional test for :py:mod`yeti.bin.findjuncs`"""
    for x in execute_helper(test_info,findjuncs_tests):
        yield x
