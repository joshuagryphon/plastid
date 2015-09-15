#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.counts_in_region`"""
import tempfile
import os

from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from plastid.util.services.decorators import catch_stderr
from plastid.test.functional.base import execute_helper
from plastid.test.ref_files import RPATH, REF_FILES, \
                                             COUNT_OPTIONS, \
                                             MASK_OPTIONS                                    
from plastid.bin.counts_in_region import main


test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "plastid.bin.counts_in_region",
    "ref_file_path"  : resource_filename("plastid","test/data/command_line"),
    "temp_file_path" : tempfile.mkdtemp(prefix="counts_in_region"),
}

# Define tests as tuples of:
# 1. Command-line style arguments to pass to ``main``
# 2. A list of reference files that output should be compared against
# 3. A list of output files created by running ``main`` with the arguments provided in (1)
# 4. A list of strings specifying how equality should be evaluated
_mask_outfile    = os.path.join(test_info["temp_file_path"],"test_counts_in_region_mask.txt")
_no_mask_outfile = os.path.join(test_info["temp_file_path"],"test_counts_in_region_no_mask.txt")
_annotation_options = "--annotation_format BED --annotation_files %s" % REF_FILES["yeast_mini_bed"]
counts_in_region_tests = [
    # test with crossmap mask
    ("%s %s %s %s" % (_mask_outfile,
                      COUNT_OPTIONS,
                      MASK_OPTIONS,
                      _annotation_options
     ),
     [REF_FILES["yeast_counts_in_region_mask"]],
     [_mask_outfile],
     ["--sort_keys region_name"]),
    # test without crossmap mask
    ("%s %s %s" % (_no_mask_outfile,
                      COUNT_OPTIONS,
                      _annotation_options
     ),
     [REF_FILES["yeast_counts_in_region_no_mask"]],
     [_no_mask_outfile],
     ["--sort_keys region_name"]),]
"""Functional tests of :py:mod:`plastid.bin.counts_in_region`.

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
    """Perform functional test for plastid.bin.counts_in_region"""
    for x in execute_helper(test_info,counts_in_region_tests):
        yield x
