#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.get_count_vectors`"""
import tempfile
import os

from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from plastid.test.functional.base import execute_helper
from plastid.test.ref_files import RPATH, REF_FILES, \
                                              COUNT_OPTIONS, \
                                              ANNOTATION_OPTIONS, \
                                              MASK_OPTIONS  
from plastid.bin.get_count_vectors import main
from plastid.util.services.decorators import catch_stderr
from plastid.util.services.mini2to3 import cStringIO


#===============================================================================
# INDEX: global constants used by tests
#===============================================================================

get_count_vectors_test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "plastid.bin.get_count_vectors",
    "ref_file_path"  : resource_filename("plastid","test/data/command_line/gen_count_vectors"),
    "temp_file_path" : tempfile.mkdtemp(prefix="get_count_vectors"),
}

# Define tests as tuples of:
# 1. Command-line style arguments to pass to ``main``
# 2. A list of reference files that output should be compared against
# 3. A list of output files created by running ``main`` with the arguments provided in (1)
# 4. A list of strings specifying how equality should be evaluated
ref_files = []
out_files = []
for line in open(REF_FILES["yeast_mini_bed"]):
    fn = line.split("\t")[3]
    ref_files.append(os.path.join(get_count_vectors_test_info["ref_file_path"], "%s.txt" % fn))            
    out_files.append(os.path.join(get_count_vectors_test_info["temp_file_path"],"%s.txt" % fn))

get_count_vectors_tests = [
    ( get_count_vectors_test_info["temp_file_path"] + COUNT_OPTIONS + MASK_OPTIONS + " --annotation_format BED --annotation_file " + REF_FILES["yeast_mini_bed"],
    ref_files,
    out_files,
    ["--no_header"]*len(ref_files)
    )
 ]
"""Functional tests of :py:mod:`plastid.bin.get_count_vectors`.

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
    for x in execute_helper(get_count_vectors_test_info,get_count_vectors_tests):
        yield x

