#!/usr/bin/env python
"""Test suite for :py:mod:`yeti.bin.metagene`"""
import tempfile
import os
from pkg_resources import resource_filename, cleanup_resources
from nose.plugins.attrib import attr
from yeti.test.functional.base import execute_helper
from yeti.test.ref_files import RPATH, REF_FILES, \
                                             COUNT_OPTIONS, \
                                             ANNOTATION_OPTIONS, \
                                             MASK_OPTIONS  
from yeti.bin.test_table_equality import main as table_test                               
from yeti.bin.metagene import main
from yeti.util.services.decorators import catch_stderr

#===============================================================================
# INDEX: global constants used by tests
#===============================================================================

TEST_INFO = { "test_method"   : catch_stderr()(main),
              "module_name"   : "yeti.bin.metagene",
              "ref_file_path"  : resource_filename("yeti","test/data/command_line"),
              "temp_file_path" : tempfile.mkdtemp(prefix="metagene"),
             }

_basename = os.path.join(TEST_INFO["temp_file_path"],"test_metagene")

#===============================================================================
# INDEX: tests
#===============================================================================

tests = [
    # test cds start
    (
        "generate %s_cds_start --downstream 100 %s %s" % (_basename, ANNOTATION_OPTIONS, MASK_OPTIONS),
        [
            REF_FILES["yeast_metagene_cds_start"],
            REF_FILES["yeast_metagene_cds_start_bed"]
        ],
        [
            _basename+"_cds_start_rois.txt",
            _basename+"_cds_start_rois.bed",
        ],
        ["","--no_header"]
    ),
    # test cds stop
    (
        "generate %s_cds_stop --upstream 100 --landmark cds_stop %s %s" % (_basename, ANNOTATION_OPTIONS, MASK_OPTIONS),
        [
            REF_FILES["yeast_metagene_cds_stop"],
            REF_FILES["yeast_metagene_cds_stop_bed"],
        ],
        [
            _basename+"_cds_stop_rois.txt",
            _basename+"_cds_stop_rois.bed",
        ],
        ["","--no_header"]
    ),
    # test count cds start
    (
        "count %s %s_cds_start --norm_region 70 150 %s" % (REF_FILES["yeast_metagene_cds_start"],
                                                    _basename,
                                                    COUNT_OPTIONS),
    [
        REF_FILES["yeast_metagene_cds_start_profile"],
        REF_FILES["yeast_metagene_cds_start_normcounts"],
        REF_FILES["yeast_metagene_cds_start_rawcounts"],
    ],
    [
        _basename+"_cds_start_metagene_profile.txt",
        _basename+"_cds_start_normcounts.txt",
        _basename+"_cds_start_rawcounts.txt"
    ],
    ["","--no_header","--no_header"]
    ),
    # test count cds stop
    (
        "count %s %s_cds_stop --norm_region 0 80 %s" % (REF_FILES["yeast_metagene_cds_stop"],
                                                     _basename,
                                                     COUNT_OPTIONS),
    [
        REF_FILES["yeast_metagene_cds_stop_profile"],
        REF_FILES["yeast_metagene_cds_stop_normcounts"],
        REF_FILES["yeast_metagene_cds_stop_rawcounts"],
    ],
    [
        _basename+"_cds_stop_metagene_profile.txt",
        _basename+"_cds_stop_normcounts.txt",
        _basename+"_cds_stop_rawcounts.txt"
    ],
    ["","--no_header","--no_header"]
    ),
]
"""Functional tests of :py:mod:`yeti.bin.metagene`.

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
    for x in execute_helper(TEST_INFO,tests):
        yield x

