#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.reformat_transcripts`"""
import tempfile
import os
import itertools
import shlex
import shutil
from nose.tools import assert_equal, assert_set_equal
from nose.plugins.attrib import attr
from pkg_resources import resource_filename, cleanup_resources
from plastid.readers.bed import BED_Reader
from plastid.readers.gff import GTF2_TranscriptAssembler
from plastid.test.functional.base import execute_helper
from plastid.test.ref_files import REF_FILES, MINI
from plastid.bin.reformat_transcripts import main
from plastid.util.services.decorators import catch_stderr

#===============================================================================
# INDEX: global constants
#===============================================================================

test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "plastid.bin.gff_to_gtf",
    "ref_file_path"  : resource_filename("plastid","test/data/command_line"),
    "temp_file_path" : tempfile.mkdtemp(prefix="reformat_transcripts"),
}


_INPUT_FORMATS  = ["BED","BigBed","GTF2","GFF3"]
_OUTPUT_FORMATS = ["BED","GTF2"]
_in_out = list(itertools.product(_INPUT_FORMATS,_OUTPUT_FORMATS))
_outfilename = os.path.join(test_info["temp_file_path"],"test_reformat_transcripts")

_test_args = { "GTF2" : " --exclude 5 8 --sort_keys 0 3 4", # exclude score and attr, which are not guaranteed to be sorted; sort on region
               "BED"  : " --sort_keys 3", # sort on region
             }

_reader_funcs = { "GTF2" : (GTF2_TranscriptAssembler, { "add_three_for_stop" : False}),
                  "BED"  : (BED_Reader,{})
                }

_default_scores = { ("BED", "BED" )   : 0.0,
                    ("BED","GTF2")    : 0.0,
                    ("GTF2","BED" )   : 0.0,
                    ("GTF2","GTF2")   : ".",
                    ("GFF3","BED")    : 0.0,
                    ("GFF3","GTF2")   : ".",
                    ("BigBed","BED")  : 0.0,
                    ("BigBed","GTF2") : 0.0,
                  }

def check_files_applying_sort(ref_file,test_file,infmt,outfmt):
    opener, args = _reader_funcs[outfmt]
    ref_transcripts  = sorted(list(opener(open(ref_file),**args)))
    test_transcripts = sorted(list(opener(open(test_file),**args)))
    
    assert_equal(len(ref_transcripts),len(test_transcripts),"%s to %s: Length mismatch in discovered transcripts. Expected '%s'. Found '%s'" % (infmt,outfmt,len(ref_transcripts),
                                                                                                                                  len(test_transcripts)))
    for tx1, tx2 in zip(ref_transcripts,test_transcripts):
        assert_equal(tx1.get_name(),tx2.get_name(),"%s to %s: Found unordered transcripts. Expected '%s'. Found '%s"'' %(infmt,outfmt,tx1.get_name(),tx2.get_name()))
        
        set1 = tx1.get_position_set()
        set2 = tx2.get_position_set()
        assert_set_equal(set1,set2,"%s to %s: Difference in position sets. Expected '%s'. Found '%s'" % (infmt,outfmt,set1,set2))

        ref_score   = str(_default_scores[(infmt,outfmt)])
        found_score = str(tx2.attr["score"])
        assert_equal(found_score,ref_score,"%s to %s: Did not find expected score. Expected: '%s'. Found '%s'" % (infmt,outfmt,ref_score,found_score))

        # BED preserves fewer attributes than GTF, so we can only test on common keys
        # We exclude "gene_id" for BEd and BigBed input because this will not match 
        # We exclude "score" because we already tested it
        #
        # by testing attr we are also implicitly testing cds_genome_end and cds_genome_start
        attr1 = tx1.attr
        attr2 = tx2.attr
        keyset = set(attr1.keys()) & set(attr2.keys()) - { "score" }
        if infmt in ("BED","BigBed") and outfmt != "BED":
            keyset -= { "gene_id" }
            
        for k in keyset:
            assert_equal(attr1[k],attr2[k],"%s to %s: Difference in attribute %s. Expected '%s'. Found '%s'" % (infmt,outfmt,k,attr1[k],attr2[k]))
        #assert_dict_equal(attr1,attr2,"%s to %s: Difference in attributes. Expected %s. Found %s" % (infmt,outfmt,attr1,attr2))



@attr(test="functional")
@attr(speed="slow")
def do_test():
    """Perform functional test for :py:mod`plastid.bin.reformat_transcripts`"""
    for infmt, outfmt in _in_out:
            
        my_outfilename = _outfilename + ("_%s_to_%s.%s" % (infmt.lower(),outfmt.lower(),outfmt.lower()))
        infile      = REF_FILES["100transcripts_%s" % infmt.lower().replace("2","").replace("3","")]
        ref_outfile = REF_FILES["100transcripts_%s" % outfmt.lower().replace("2","").replace("3","")]
        argstr = "%s --annotation_files %s --annotation_format %s --output_format %s" % (my_outfilename,infile,infmt,outfmt)
        test_info["test_method"](shlex.split(argstr))
        yield check_files_applying_sort, ref_outfile, my_outfilename, infmt, outfmt


def tearDownModule():
    if test_info["temp_file_path"] != "":
        shutil.rmtree(test_info["temp_file_path"])
