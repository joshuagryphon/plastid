#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.findjuncs`"""
import tempfile
import os

from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from plastid.test.functional.base import execute_helper
from plastid.test.ref_files import REF_FILES
from plastid.bin.findjuncs import main
from plastid.util.services.decorators import catch_stderr


#===============================================================================
# INDEX: global constants used by tests
#===============================================================================

test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "plastid.bin.findjuncs",
    "ref_file_path"  : resource_filename("plastid","test/data/"),
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
"""Functional tests of :py:mod:`plastid.bin.findjuncs`.

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
    """Perform functional test for :py:mod`plastid.bin.findjuncs`"""
    for x in execute_helper(test_info,findjuncs_tests):
        yield x
