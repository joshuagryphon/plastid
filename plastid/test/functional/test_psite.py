#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.psite`"""
import tempfile
import os

from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from plastid.util.services.decorators import catch_stderr
from plastid.test.functional.base import execute_helper
from plastid.test.ref_files import REF_FILES
from plastid.bin.psite import main

#===============================================================================
# INDEX: global constants
#===============================================================================

test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "plastid.bin.psite",
    "ref_file_path"  : resource_filename("plastid","test/data/command_line"),
    "temp_file_path" : tempfile.mkdtemp(prefix="psite"),
}

_outbase = os.path.join(test_info["temp_file_path"],"test_psite")

#===============================================================================
# INDEX: tests
#===============================================================================

psite_tests = [
    ("%s %s --count_files %s --norm_region 70 150 --require_upstream --min_length 26 --max_length 31"  % (REF_FILES["yeast_metagene_cds_start"],
                                                                                                          _outbase,
                                                                                                          REF_FILES["yeast_rp_bam"],)  ,
     [REF_FILES["yeast_psite"]],
     [_outbase+"_p_offsets.txt"],
     [""]
    ),
    ("%s %s --count_files %s --norm_region 70 150 --require_upstream --min_length 26 --max_length 31 --constrain 0 11"  % (REF_FILES["yeast_metagene_cds_start"],
                                                                                                          _outbase,
                                                                                                          REF_FILES["yeast_rp_bam"],)  ,
     [REF_FILES["yeast_psite_constrain"]],
     [_outbase+"_p_offsets.txt"],
     [""]
    ),
    # test using --normalize_over instead of --norm_region
    ("%s %s --count_files %s --normalize_over 20 100 --require_upstream --min_length 26 --max_length 31"  % (REF_FILES["yeast_metagene_cds_start"],
                                                                                                          _outbase,
                                                                                                          REF_FILES["yeast_rp_bam"],)  ,
     [REF_FILES["yeast_psite"]],
     [_outbase+"_p_offsets.txt"],
     [""]
    ),
    # test using --normalize_over instead of --norm_region
    ("%s %s --count_files %s --normalize_over 20 100 --require_upstream --min_length 26 --max_length 31 --constrain 0 11"  % (REF_FILES["yeast_metagene_cds_start"],
                                                                                                          _outbase,
                                                                                                          REF_FILES["yeast_rp_bam"],)  ,
     [REF_FILES["yeast_psite_constrain"]],
     [_outbase+"_p_offsets.txt"],
     [""]
    ),               
]


# /home/joshua/projects/plastid/plastid/test/data/command_line/gen_cds_start_rois.txt /tmp/psite5SY5Wc/test_psite --count_files /home/joshua/projects/plastid/plastid/test/data/command_line/gen_reads.bam --norm_region 70 150 --require_upstream --min_length 26 --max_length 31

"""Functional tests of :py:mod:`plastid.bin.psite`.

Tests are specified as tuples of:

    1. Command-line style arguments to pass to :py:func:`main`

    2. A list of reference files that output should be compared against

    3. A list of output files created by running :py:func:`main`
       with the arguments provided in (1)

    4. A list of strings specifying how equality should be evaluated
"""

#===============================================================================
# INDEX: test/helper functions
#===============================================================================

@attr(test="functional")
@attr(speed="slow")
def do_test():
    """Perform functional test for :py:mod:`plastid.bin.psite`"""
    for x in execute_helper(test_info,psite_tests):
        yield x