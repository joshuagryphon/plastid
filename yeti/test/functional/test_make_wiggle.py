#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.make_wiggle`"""
import tempfile
import shutil
import os
import shlex
import nose

from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from plastid.test.functional.base import execute_helper
from plastid.test.ref_files import RPATH, REF_FILES, \
                                             COUNT_OPTIONS, \
                                             ANNOTATION_OPTIONS, \
                                             MASK_OPTIONS  
from plastid.genomics.genome_array import GenomeArray                                             
from plastid.bin.make_wiggle import main
from plastid.util.services.decorators import catch_stderr

#===============================================================================
# INDEX: global constants used by tests
#===============================================================================

test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "plastid.bin.make_wiggle",
    "ref_file_path"  : resource_filename("plastid","test/data/command_line"),
    "temp_file_path" : tempfile.mkdtemp(prefix="make_wiggle"),
}
"""Constants used in tests below"""

#===============================================================================
# INDEX: programmatically define tests
#===============================================================================

# Define tests as tuples of:
# 1. Command-line style arguments to pass to ``main``
# 2. A list of reference files that output should be compared against
# 3. A list of output files created by running ``main`` with the arguments provided in (1)
# 4. A list of strings specifying how equality should be evaluated
make_wiggle_tests = []
samples = { "fiveprime_variable" : "--fiveprime_variable --offset %s" % REF_FILES["yeast_psite"],
          }

for side in ("fiveprime","threeprime"):
    for offset in (0,15):
        samples["%s_%s" % (side,offset)] = "--%s --offset %s" % (side,offset)

for nibble in (0,12):
    samples["center_%s" % nibble] = "--center --nibble %s" % nibble

for k,v in samples.items():
    fbase = os.path.join(test_info["temp_file_path"],"test_make_wiggle_%s" % k)
    make_wiggle_tests.append( ("%s --count_files %s -o %s " % (v,REF_FILES["yeast_rp_bam"],fbase),
                       [os.path.join(test_info["ref_file_path"],"gen_reads_%s_fw.wig" % k),
                        os.path.join(test_info["ref_file_path"],"gen_reads_%s_rc.wig" % k),
                       ],
                       [os.path.join(test_info["temp_file_path"],"test_make_wiggle_%s_fw.wig" % k),
                        os.path.join(test_info["temp_file_path"],"test_make_wiggle_%s_rc.wig" % k),
                       ]
                       )
                     )


#===============================================================================
# INDEX: test functions and helpers 
#===============================================================================

def do_check(argstr,ref_files,test_files):
    """Check that test wiggle files match reference files
    
    Parameters
    ----------
    ref_files : tuple
        Tuple of filenames corresponding to forward and reverse strand wiggles
        for reference dataset
    
    test_files : list
        Tuple of filenames corresponding to forward and reverse strand wiggles
        for test dataset
    """
    err_message = "Unequal output in files %s vs %s, module %s" % (ref_files,test_files,test_info["module_name"])
    test_info["test_method"](shlex.split(argstr))
    ref_ga  = GenomeArray()
    test_ga = GenomeArray()
    ref_ga.add_from_wiggle(open(ref_files[0]),"+")
    ref_ga.add_from_wiggle(open(ref_files[1]),"-")
    test_ga.add_from_wiggle(open(test_files[0]),"+")
    test_ga.add_from_wiggle(open(test_files[1]),"-")
    nose.tools.assert_true(ref_ga.__eq__(test_ga,tol=1e-8),err_message)

@attr(test="functional")
@attr(speed="slow")
def do_test():
    """Execute functional test for :py:mod:`plastid.bin.make_wiggle`"""
    for argstr, ref_files, test_files in make_wiggle_tests:
        yield do_check, argstr, ref_files, test_files

def tearDownModule():
    """Remove test dataset files after unit tests are complete"""
    if test_info["temp_file_path"] != "":
        shutil.rmtree(test_info["temp_file_path"])
            
    cleanup_resources()

