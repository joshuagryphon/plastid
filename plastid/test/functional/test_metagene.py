#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.metagene`"""
import tempfile
import os
from pkg_resources import resource_filename, cleanup_resources
from nose.plugins.attrib import attr
from plastid.test.functional.base import execute_helper
from plastid.test.ref_files import RPATH, REF_FILES, \
                                             COUNT_OPTIONS, \
                                             ANNOTATION_OPTIONS, \
                                             MASK_OPTIONS  
from plastid.bin.test_table_equality import main as table_test                               
from plastid.bin.metagene import main
from plastid.util.services.decorators import catch_stderr

#===============================================================================
# INDEX: global constants used by tests
#===============================================================================

TEST_INFO = { "test_method"    : catch_stderr()(main),
              "module_name"    : "plastid.bin.metagene",
              "ref_file_path"  : resource_filename("plastid","test/data/command_line"),
              "temp_file_path" : tempfile.mkdtemp(prefix="metagene"),
             }

_basename = os.path.join(TEST_INFO["temp_file_path"],"test_metagene")

#===============================================================================
# INDEX: tests
#===============================================================================

tests = [
    # test generate cds start
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
    # test generate cds stop
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
    # test count cds start with --norm_region
    (
        "count %s %s_cds_start --keep --norm_region 70 150 %s" % (REF_FILES["yeast_metagene_cds_start"],
                                                    _basename,
                                                    COUNT_OPTIONS),
    [
        REF_FILES["yeast_metagene_cds_start_profile"],
        REF_FILES["yeast_metagene_cds_start_normcounts"],
        REF_FILES["yeast_metagene_cds_start_rawcounts"],
    ],
    [
        _basename+"_cds_start_metagene_profile.txt",
        _basename+"_cds_start_normcounts.txt.gz",
        _basename+"_cds_start_rawcounts.txt.gz"
    ],
    ["","--no_header","--no_header"]
    ),
    # test count cds stop with --norm_region
    (
        "count %s %s_cds_stop --keep --norm_region 0 80 %s" % (REF_FILES["yeast_metagene_cds_stop"],
                                                     _basename,
                                                     COUNT_OPTIONS),
    [
        REF_FILES["yeast_metagene_cds_stop_profile"],
        REF_FILES["yeast_metagene_cds_stop_normcounts"],
        REF_FILES["yeast_metagene_cds_stop_rawcounts"],
    ],
    [
        _basename+"_cds_stop_metagene_profile.txt",
        _basename+"_cds_stop_normcounts.txt.gz",
        _basename+"_cds_stop_rawcounts.txt.gz"
    ],
    ["","--no_header","--no_header"]
    ),
    # test count cds start, using --normalize_over
    (
        "count %s %s_cds_start --keep --normalize_over 20 100 %s" % (REF_FILES["yeast_metagene_cds_start"],
                                                    _basename,
                                                    COUNT_OPTIONS),
    [
        REF_FILES["yeast_metagene_cds_start_profile"],
        REF_FILES["yeast_metagene_cds_start_normcounts"],
        REF_FILES["yeast_metagene_cds_start_rawcounts"],
    ],
    [
        _basename+"_cds_start_metagene_profile.txt",
        _basename+"_cds_start_normcounts.txt.gz",
        _basename+"_cds_start_rawcounts.txt.gz"
    ],
    ["","--no_header","--no_header"]
    ),
    # test count cds stop, using --normalize_over
    (
        "count %s %s_cds_stop --keep --normalize_over '-100' '-20' %s" % (REF_FILES["yeast_metagene_cds_stop"],
                                                     _basename,
                                                     COUNT_OPTIONS),
    [
        REF_FILES["yeast_metagene_cds_stop_profile"],
        REF_FILES["yeast_metagene_cds_stop_normcounts"],
        REF_FILES["yeast_metagene_cds_stop_rawcounts"],
    ],
    [
        _basename+"_cds_stop_metagene_profile.txt",
        _basename+"_cds_stop_normcounts.txt.gz",
        _basename+"_cds_stop_rawcounts.txt.gz"
    ],
    ["","--no_header","--no_header"]
    ),         
]
"""Functional tests of :py:mod:`plastid.bin.metagene`.

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

