#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.cs`"""
import tempfile
import os

from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from plastid.test.functional.base import execute_helper
from plastid.test.ref_files import COUNT_OPTIONS, \
                                   ANNOTATION_OPTIONS, \
                                   MASK_OPTIONS  
from plastid.bin.cs import main

from plastid.util.services.decorators import catch_stderr


#===============================================================================
# INDEX: global constants used by tests
#===============================================================================

test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "plastid.bin.cs",
    "ref_file_path"  : resource_filename("plastid","test/data/command_line"),
    "temp_file_path" : tempfile.mkdtemp(prefix="cs"),
}
"""Constants used by multiple tests"""

_outbase = os.path.join(test_info["temp_file_path"],"test_cs")
_file_stubs = [
                'merged.txt',
                'gene.positions',
                'transcript.positions',

                'gene_cds.bed',
                'gene_exon.bed',
                'gene_masked.bed',
                'gene_utr3.bed',
                'gene_utr5.bed',
                'transcript_cds.bed',
                'transcript_exon.bed',
                'transcript_masked.bed',
                'transcript_utr3.bed',
                'transcript_utr5.bed',
              ]


#===============================================================================
# INDEX: tests to execute
#===============================================================================

cs_tests = [
    ("generate %s %s %s" % (_outbase,ANNOTATION_OPTIONS,MASK_OPTIONS),
     [os.path.join(test_info["ref_file_path"],  "gen_cs_"+X) for X in _file_stubs],
     [os.path.join(test_info["temp_file_path"],"test_cs_"+X) for X in _file_stubs],
     ["--no_header","--sort_keys region","--sort_keys region",] + ["--no_header --sort_keys 3"]*10
    ),
    ("count %s_gene.positions %s_count %s " % (_outbase,_outbase,COUNT_OPTIONS),
     [os.path.join(test_info["ref_file_path"],'gen_cs_count_gene.txt')],
     [os.path.join(test_info["temp_file_path"],"test_cs_count.txt")],
     ["--sort_keys region"]
    ),
    ("count %s_gene.positions %s_count %s --sum 1e9" % (_outbase,_outbase,COUNT_OPTIONS),
     [os.path.join(test_info["ref_file_path"],'gen_cs_count_gene_sum_1.txt')],
     [os.path.join(test_info["temp_file_path"],"test_cs_count.txt")],
     ["--sort_keys region"]
    ),
    #("",[],[],[]),
]
"""Functional tests of :py:mod:`plastid.bin.cs`.

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
    """Perform functional test for :py:mod:`plastid.bin.cs`"""
    for x in execute_helper(test_info,cs_tests):
        yield x
