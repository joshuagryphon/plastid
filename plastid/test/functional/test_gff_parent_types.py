#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.gff_parent_types`"""
import tempfile
import os

from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from plastid.test.functional.base import execute_helper
from plastid.test.ref_files import REF_FILES
from plastid.bin.gff_parent_types import main
from plastid.util.services.decorators import catch_stderr


#===============================================================================
# INDEX: global constants used by tests
#===============================================================================

test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "plastid.bin.gff_parent_types",
    "ref_file_path"  : resource_filename("plastid","test/data/command_line"),
    "temp_file_path" : tempfile.mkdtemp(prefix="test_gff_parent_types"),
}

# Define tests as tuples of:
# 1. Command-line style arguments to pass to ``main``
# 2. A list of reference files that output should be compared against
# 3. A list of output files created by running ``main`` with the arguments provided in (1)
# 4. A list of strings specifying how equality should be evaluated
_outfile = os.path.join(test_info["temp_file_path"],"test_parent_types.txt")
gff_parent_types_tests = [
                         # test ot ignoring any feature types 
                         ("%s %s" % (REF_FILES["yeast_gff"],_outfile),
                          [REF_FILES["yeast_parent_child"]],
                          [_outfile],
                          [""]),
                         # test ignoring manu feature types
                         ("%s %s --exclude ARS tRNA LTR_retrotransposon StopFeature X_element" % (REF_FILES["yeast_gff"],_outfile),
                          [REF_FILES["yeast_parent_child_exclude_cols"]],
                          [_outfile],
                          [""]),
]

@attr(test="functional")
@attr(speed="slow")
def do_test():
    """Perform functional test for :py:mod`plastid.bin.gff_parent_types`"""
    for x in execute_helper(test_info,gff_parent_types_tests):
        yield x
