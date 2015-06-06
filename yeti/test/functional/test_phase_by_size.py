#!/usr/bin/env python
"""Test suite for :py:mod:`yeti.bin.phase_by_size`"""
import tempfile
import os

from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from yeti.util.services.decorators import catch_stderr
from yeti.test.functional.base import execute_helper
from yeti.test.ref_files import RPATH, REF_FILES, \
                                             COUNT_OPTIONS, \
                                             ANNOTATION_OPTIONS, \
                                             MASK_OPTIONS  
from yeti.bin.phase_by_size import main
from yeti.util.services.mini2to3 import cStringIO

#===============================================================================
# INDEX: global constants
#===============================================================================

test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "yeti.bin.phase_by_size",
    "ref_file_path"  : resource_filename("yeti","test/data/command_line"),
    "temp_file_path" : tempfile.mkdtemp(prefix="phase_by_size"),
}

# Define tests as tuples of:
# 1. Command-line style arguments to pass to ``main``
# 2. A list of reference files that output should be compared against
# 3. A list of output files created by running ``main`` with the arguments provided in (1)
# 4. A list of strings specifying how equality should be evaluated
_basename = os.path.join(test_info["temp_file_path"],"test_phase")

#===============================================================================
# INDEX: tests
#===============================================================================

phase_by_size_tests = [
    ("%s --min_length 26 --max_length 31 --annotation_files %s --annotation_format BED %s" % (_basename,
                                                                                             REF_FILES["yeast_mini_bed"],
                                                                                             COUNT_OPTIONS),
     [REF_FILES["yeast_phasing"]],
     [_basename+"_phasing.txt"],
     [""]
    ),
]
"""Functional tests of :py:mod:`yeti.bin.phase_by_size`.

Tests are specified as tuples of:

    1. Command-line style arguments to pass to :py:func:`main`

    2. A list of reference files that output should be compared against

    3. A list of output files created by running :py:func:`main`
       with the arguments provided in (1)

    4. A list of strings specifying how equality should be evaluated
"""


#===============================================================================
# INDEX: test functions
#===============================================================================

@attr(test="functional")
@attr(speed="slow")
def do_test():
    """Perform functional test for :py:mod:`yeti.bin.phase_by_size`"""
    for x in execute_helper(test_info,phase_by_size_tests):
        yield x

