#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.slidejuncs`"""
import tempfile
import shutil
import os
import itertools
import numpy
import copy
import shlex

from nose.tools import assert_set_equal
from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from plastid.test.ref_files import REF_FILES
from plastid.util.io.filters import CommentReader
from plastid.readers.bed import BED_to_SegmentChain
from plastid.bin.slidejuncs import main
from plastid.util.services.decorators import catch_stderr, catch_stdout

#===============================================================================
# INDEX: global constants
#===============================================================================

test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "plastid.bin.slidejuncs",
    "ref_file_path"  : resource_filename("plastid","test/data/command_line"),
    "temp_file_path" : tempfile.mkdtemp(prefix="slidejuncs"),
}

OUTBASE = os.path.join(test_info["temp_file_path"],"slidejuncs_%s_%s_%s")

# reference files used
FASTA_FILE   = REF_FILES["slidejuncs_seqs"]
COMBINED_BED = REF_FILES["slidejuncs_input"]
REF_FILE     = REF_FILES["slidejuncs_ref"]
MASK_FILE    = REF_FILES["slidejuncs_crossmap"]




#===============================================================================
# INDEX: helper functions
#===============================================================================

def get_junction_names_from_file(filename):
    """Return a set of names of splice junctions in a file
    
    Parameters
    ----------
    filename : str
        Fully-qualified path to file
    
    Returns
    -------
    set
        Set of splice junction names found in ``filename``    
    """
    return set([X.spanning_segment.chrom for X in BED_to_SegmentChain(CommentReader(open(filename)))])
    
def get_junction_names_from_category(cat):
    """Helper function to deduce names of splice junctions/genes in each
    junction category from input files
    
    Parameters
    ----------
    cat : str
        Category of splice junction. Must be present in ``junction_categories``
    
    Returns
    -------
    set
        Set of junction names in input splice junction category ``cat``
    """
    return get_junction_names_from_file(REF_FILES["slidejuncs_%s" % cat])


#===============================================================================
# INDEX: programmatic definition of tests and expected results
#===============================================================================

# query junctions
junction_categories = ["known_juncs_non_crossmap",
                       "known_juncs_crossmap",
                       "to_slide_known_non_crossmap",
                       "to_slide_known_crossmap",
                       "noncan_no_ref",
                       "expected_untouched",
                       ]
"""Categories of splice junctions taken as input in tests below"""

junctions = { K : get_junction_names_from_category(K) for K in junction_categories }
"""Dictionary mapping categories of splice junctions to names of junctions
used in tests below"""

# classification order is repetitive > reference > canonical > untouched

# set up tests as configurations of kwargs
# we'll represent these as 3-tuples, where 1 in a position indicates whether
# or not we'll supply:
#     0. a crossmap
#     1. a reference file
#     2. the --slide_canonical flag
#     3. always value 1 (untouched state)
test_options = list(itertools.product((0,1),repeat=3))
test_options = [numpy.array(list(X) + [1]) for X in test_options]
base_command_line = " ".join([COMBINED_BED,OUTBASE,"--sequence_file %s" % FASTA_FILE])
switches = [" --mask_annotation_format BED --mask_annotation_file %s" % MASK_FILE,
            " --ref %s" % REF_FILE,
            " --slide_canonical",
            ""]
"""Command-line switches used for building tests"""

labels = ["repetitive",
          "shifted_known",
          "shifted_canonical",
          "untouched"]
"""Labels for various junctions under different test conditions"""
# Classification preferences for each set of query junctions
#
# 1 represents the ability to be classified in a category denoted by position,
#
#     0. repetitive
#     1. reference
#     2. canonical
#     3. untouched
#
# These positions correspond to the test_option positions above in the sense
# that if both test_option[position] and classification_pref[k][position] are 1,
# then junctions in group k can be classified with the label corresponding
# to that position. In practice, the leftmost position in which both bits are 1
# should be chosen by the script
classification_prefs = { "known_juncs_non_crossmap"    : (0,1,1,1),
                         "known_juncs_crossmap"        : (1,1,1,1),
                         "to_slide_known_non_crossmap" : (0,1,1,1),
                         "to_slide_known_crossmap"     : (1,1,1,1),
                         "noncan_no_ref"               : (0,0,1,1),
                         "expected_untouched"          : (0,0,0,1),
                    }
classification_prefs = { K : numpy.array(V) for K,V in classification_prefs.items() }
"""Classification preferences for each category of query junction"""

output_files = { K : OUTBASE+("_%s.bed" % K) for K in labels }
"""Output files written by each test"""

# tests will be defined as tuples of command-line arguments,
# and dictionaries mapping classifications to splice junctions that
# should fall into that category for a given run
tests = []
for my_opts in test_options:
    # build command-line arguments
    stmp = base_command_line + " ".join([switches[X] for X in my_opts.nonzero()[0]])
    stmp = stmp % (my_opts[0],my_opts[1],my_opts[2])
    result_dict = {}
    
    # determine results
    for k,v in classification_prefs.items():
        label = labels[(v & my_opts).argmax()]
        try:
            result_dict[label] |= junctions[k]
        except KeyError:
            result_dict[label] = copy.deepcopy(junctions[k])
    
    tests.append((my_opts,stmp,result_dict))


#===============================================================================
# INDEX: tests
#===============================================================================

def compare_sets(label,found,expected,crossmap,ref,slide):
    message = "Failed test %s.\n\tCrossmap: %s\n\tRef: %s\n\tslide: %s\nItems expected but not found: %s\nVice versa: %s\n" % (label,
                                                                                                                            crossmap,
                                                                                                                            ref,
                                                                                                                            slide,
                                                                                                                            expected-found,
                                                                                                                            found-expected)
    assert_set_equal(expected,found,message)
    
@attr(test="functional")
def do_test():
    for tup, argstr, expected_results in tests:
        test_info["test_method"](shlex.split(argstr))
        
        for label, filename in output_files.items():
            found    = get_junction_names_from_file(filename % (tup[0], tup[1], tup[2]))
            expected = expected_results.get(label,set([]))
            yield compare_sets, label, found, expected, tup[0], tup[1], tup[2]
 
def tearDownModule():
    """Remove test dataset files after unit tests are complete"""
    if test_info["temp_file_path"] != "":
        shutil.rmtree(test_info["temp_file_path"])
             
    cleanup_resources()

