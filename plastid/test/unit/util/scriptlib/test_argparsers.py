#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.util.scriptlib.argparsers`"""
import unittest
import shlex
import argparse
import copy
import numpy
import itertools
import warnings
from pkg_resources import cleanup_resources
from twobitreader import TwoBitFile
from nose.plugins.attrib import attr
from nose.tools import assert_equal, assert_greater, assert_raises,  \
    assert_in, assert_not_in, assert_true, assert_list_equal

from plastid.test.ref_files import MINI, REF_FILES
from plastid.util.io.filters import CommentReader
from plastid.readers.bed import BED_to_SegmentChain
from plastid.util.services.exceptions import MalformedFileError
from plastid.genomics.roitools import SegmentChain, Transcript
from plastid.genomics.map_factories import FivePrimeMapFactory,\
    ThreePrimeMapFactory, CenterMapFactory, VariableFivePrimeMapFactory
from plastid.genomics.genome_array import GenomeArray, BAMGenomeArray, \
    SparseGenomeArray
from plastid.genomics.genome_hash import GenomeHash, BigBedGenomeHash

from plastid.util.scriptlib.argparsers import PrefixNamespaceWrapper,\
                                           get_alignment_file_parser,\
                                           get_genome_array_from_args,\
                                           get_annotation_file_parser,\
                                           get_transcripts_from_args,\
                                           get_segmentchain_file_parser,\
                                           get_segmentchains_from_args,\
                                           get_mask_file_parser,\
                                           get_genome_hash_from_mask_args,\
                                           get_sequence_file_parser,\
                                           get_seqdict_from_args,\
                                           _parse_variable_offset_file

from plastid.util.services.mini2to3 import StringIO

warnings.simplefilter("ignore")


def tearDownModule():
    cleanup_resources()

@attr(test="unit")
class TestParseVariableOffsetFile(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.ref_file = REF_FILES["yeast_psite"]
        cls.lines = CommentReader(open(cls.ref_file)).read().strip("\n").split("\n")
        cls.expected = { 26 : 12,
                         27 : 12,
                         28 : 13,
                         29 : 13,
                         30 : 14,
                         31 : 13,
                         "default" : 13,
                        }
    
    def test_no_default(self):
        my_expected = copy.deepcopy(self.expected)
        my_expected.pop("default")
        my_lines = "\n".join(self.lines[:-1]) 
        dtmp = _parse_variable_offset_file(StringIO.StringIO(my_lines))
        self.assertDictEqual(my_expected,dtmp)

    def test_with_default(self):
        dtmp = _parse_variable_offset_file(CommentReader(open(self.ref_file)))
        self.assertDictEqual(self.expected,dtmp)

    def test_only_default(self):
        dtmp =  _parse_variable_offset_file(StringIO.StringIO(self.lines[-1]))
        my_expected = { "default" : 13 }
        self.assertDictEqual(my_expected,dtmp)
    
    def test_redefined_raises_exception(self):
        # make sure files with multiple entries for the same key error out
        my_lines = "\n".join(self.lines + self.lines)
        fakefile = StringIO.StringIO(my_lines)
        fakefile.__name__ = "Fake file" 
        self.assertRaises(MalformedFileError,_parse_variable_offset_file,fakefile)

    def test_wrong_cols_raises_exception(self):
        too_many = """25    12
26    13
27    14
28    15    52
""".replace("    ","\t")
        fakefile = StringIO.StringIO(too_many)
        fakefile.__name__ = "Fake file" 
        self.assertRaises(MalformedFileError,_parse_variable_offset_file,fakefile)

        too_few = """25    12
26    13
27
28    15
""".replace("    ","\t")
        fakefile = StringIO.StringIO(too_few)
        fakefile.__name__ = "Fake file" 
        self.assertRaises(MalformedFileError,_parse_variable_offset_file,fakefile)
    
    def test_non_int_raises_exception(self):
        non_int_key = """25    12
26    13
a    14
28    15
""".replace("    ","\t")
        fakefile = StringIO.StringIO(non_int_key)
        fakefile.__name__ = "Fake file" 
        self.assertRaises(MalformedFileError,_parse_variable_offset_file,fakefile)

        non_int_val = """25    12
26    13
27    q
28    15
""".replace("    ","\t")
        fakefile = StringIO.StringIO(non_int_val)
        fakefile.__name__ = "Fake file" 
        self.assertRaises(MalformedFileError,_parse_variable_offset_file,fakefile)


        null_val = """25    12
26    13
27    
28    15
""".replace("    ","\t")
        fakefile = StringIO.StringIO(null_val)
        fakefile.__name__ = "Fake file" 
        self.assertRaises(MalformedFileError,_parse_variable_offset_file,fakefile)

 
#=============================================================================
# INDEX: constants and helper functions used for tests in multiple sections
#        below.
#=============================================================================

offset_file    = REF_FILES["yeast_psite"]
bamfilename    = MINI["bamfile"]
bowtiefilename = MINI["bowtie_file"]
wigglefilename = MINI["wig_bedgraph"]
regions = list(BED_to_SegmentChain(open(MINI["bed_file"])))
spliced_regions = [X for X in regions if len(X) > 1] # set(filter(lambda x: len(x) > 1, regions))
unique_regions  = [X for X in regions if "repeat" not in X.get_name()] #set(filter(lambda x: "repeat" not in x.get_name(), regions))
introns = [X for X in regions if "intron" not in X.get_name()] #set(filter(lambda x: "intron" in x.get_name(),regions))

unspliced_unique_regions = [X for X in unique_regions if all([len(X) == 1,"intron" not in X.get_name()])] #unique_regions - spliced_regions - introns

assert_greater(len(unspliced_unique_regions),0)

alignment_file_parser_disableable = ("fiveprime",
                                "threeprime",
                                "center",
                                "offset",
                                "nibble",
                                "min_length",
                                "max_length",
                                "normalize",
                                "big_genome")
alignment_file_parser_opts = ("big_genome",
                              "count_files",
                              "countfile_format",
                              "mapping",
                              "max_length",
                              "min_length",
                              "nibble",
                              "normalize",
                              "offset")


annotation_file_parser_opts = {
    "annotation_format",
    "add_three",
    "annotation_files",
    "tabix"
}
annotation_file_gff_opts = { "gff_exon_types", "gff_transcript_types"}
annotation_file_parser_disableable = annotation_file_parser_opts - set(["annotation_files"])


ivcollection_file_parser_opts        = annotation_file_parser_opts - set(["add_three"])
ivcollection_file_parser_disableable = annotation_file_parser_disableable - set(["add_three"])

mask_file_parser_opts        = ivcollection_file_parser_opts
mask_file_parser_disableable = ivcollection_file_parser_disableable

alignment_file_parser    = get_alignment_file_parser()
annotation_file_parser   = get_annotation_file_parser()
ivcollection_file_parser = get_segmentchain_file_parser()
sequence_file_parser     = get_sequence_file_parser()
mask_file_parser         = get_mask_file_parser()

def check_prefix(parser_fn,opts):
    """Helper function to test prefix appending to various parsers

    Parameters
    ----------
    parser_fn : function
        Function that returns a :py:class:`~argparse.ArgumentParser`

    opts : list
        list of strings of options, without prefixes, included
        in the parser
    """
    parser = parser_fn(prefix="test")
    args = parser.parse_args([])
    for opt in opts:
        # getattr raises an AttributeError if attribute not present
        # so it is a pseudotest
        getattr(args,"%s%s" % ("test",opt))

def check_disabled_options(parser_fn,disable_opts):
    """Helper function to test disabling of options in various parsers

    Parameters
    ----------
    parser_fn : function
        Function that returns a :py:class:`~argparse.ArgumentParser`

    disable_opts : list
        list of strings of options, without prefixes, that can be 
        disabled in the parser
    """
    # kind of a hack. make sure options are disabled
    # by asserting that they don't appear in --help
    for d in disable_opts:
        parser = parser_fn(disabled=[d])
        assert_not_in("\n  --%s " % d,parser.format_help())

@attr(test="unit")
def test_parser_prefixes():
    tuples = [(get_alignment_file_parser,alignment_file_parser_opts),
              (get_annotation_file_parser,annotation_file_parser_opts),
              (get_segmentchain_file_parser,ivcollection_file_parser_opts),
              (get_mask_file_parser,mask_file_parser_opts),
              ]
    for fn, opts in tuples:
        yield check_prefix, fn, opts

@attr(test="unit")
def test_parser_disabled_options():
    tuples = [(get_alignment_file_parser,alignment_file_parser_disableable),
              (get_annotation_file_parser,annotation_file_parser_disableable),
              (get_segmentchain_file_parser,ivcollection_file_parser_disableable),
              (get_mask_file_parser,mask_file_parser_disableable),
              ]
    for fn, opts in tuples:
        yield check_disabled_options, fn, opts


#=============================================================================
# INDEX: tests for alignment file parsing
#=============================================================================


@attr(test="unit")
def test_alignment_parser_no_count_file_causes_exit():
    args = alignment_file_parser.parse_args([])
    assert_raises(SystemExit,get_genome_array_from_args,args)
    
@attr(test="unit")
def test_alignment_parser_no_mapping_causes_exit():
    argstr = "--countfile_format BAM --count_files %s" % bamfilename
    args = alignment_file_parser.parse_args(shlex.split(argstr))
    assert_raises(SystemExit,get_genome_array_from_args,args)
    
    argstr = "--countfile_format bowtie --count_files %s" % bamfilename
    args = alignment_file_parser.parse_args(shlex.split(argstr))
    assert_raises(SystemExit,get_genome_array_from_args,args)
    
@attr(test="unit")
def test_alignment_parser_open_bam_makes_bamgenomearray():
    argstr = "--countfile_format BAM --count_files %s --fiveprime" % bamfilename
    args = alignment_file_parser.parse_args(shlex.split(argstr))
     
    assert_equal(bamfilename,args.count_files[0])
     
    ga = get_genome_array_from_args(args)
    assert_true(isinstance(ga,BAMGenomeArray))
 
@attr(test="unit")
def test_alignment_parser_open_bowtie_makes_genomearray():
    argstr = "--countfile_format bowtie --count_files %s --fiveprime" % bowtiefilename
    args = alignment_file_parser.parse_args(shlex.split(argstr))
    ga = get_genome_array_from_args(args)
    assert_true(isinstance(ga,GenomeArray))
  
@attr(test="unit")
def test_alignment_parser_open_wiggle_makes_genomearray():
    argstr = "--countfile_format wiggle --count_files %s" % wigglefilename
    args = alignment_file_parser.parse_args(shlex.split(argstr))
    ga = get_genome_array_from_args(args)
    assert_true(isinstance(ga,GenomeArray))
  
@attr(test="unit")
def test_alignment_parser_bowtie_with_sparsegenomearray():
    argstr = "--countfile_format bowtie --count_files %s --fiveprime --big_genome" % bowtiefilename
    args = alignment_file_parser.parse_args(shlex.split(argstr))
    ga = get_genome_array_from_args(args)
    assert_true(isinstance(ga,SparseGenomeArray))

@attr(test="unit")
def test_alignment_parser_wiggle_with_sparsegenomearray():
    argstr = "--countfile_format wiggle --count_files %s --big_genome" % wigglefilename
    args = alignment_file_parser.parse_args(shlex.split(argstr))
     
    ga = get_genome_array_from_args(args)
    assert_true(isinstance(ga,SparseGenomeArray))

def check_genome_array_against_vector(ga,mapping,offset_or_nibble):
    """Compare values in a |GenomeArray| against known vectors
    across a set of unspliced, uniquely mapping regions

    Parameters
    ----------
    ga : |GenomeArray|, |SparseGenomeArray|, or |BAMGenomeArray|
        Genome array to test

    mapping : str
        String specifying mapping rule. "fiveprime", "threeprime",
        "center", or "fiveprime_variable"

    offset_or_nibble : int or str
        Integer specifying offset or nibble, or str specifying filename
        of variable offset dict
    """
    if mapping != "fiveprime_variable":
        vec_fw = numpy.loadtxt(MINI["vec_%s_%s_fw.txt" % (mapping,offset_or_nibble)])
        vec_rc = numpy.loadtxt(MINI["vec_%s_%s_rc.txt" % (mapping,offset_or_nibble)])
    else:
        vec_fw = numpy.loadtxt(MINI["vec_threeprime_15_fw.txt"])
        vec_rc = numpy.loadtxt(MINI["vec_threeprime_15_rc.txt"])

    for region in unspliced_unique_regions:
        vec = []
        if region.strand == "+":
            for iv in region:
                vec.extend(vec_fw[iv.start:iv.end])
        else:
            for iv in region:
                vec.extend(vec_rc[iv.start:iv.end])
            vec = vec[::-1]

        vec  = numpy.array(vec)
        vec2 = numpy.array(region.get_counts(ga))
        diff = abs(vec2 - vec)
        assert_true((diff <= 1e-8).all(),"Maximum difference exceeded on %s: %s\n%s\n%s" % (region.get_name(),diff.max(),vec[:40],vec2[:40]))

def check_bam_mapping(mapping,offset_or_nibble):
    """Generator function to test import of BAM files under many read mapping rules

    Parameters
    ----------
    mapping : str
        String specifying mapping rule. "fiveprime",
        "threeprime", "center", or "fiveprime_variable"

    offset_or_nibble : int or str
        Integer specifying offset or nibble (only for center mapping)
        or str specifying filename of variable offset dict text file
    """
    argstr = "--countfile_format BAM --count_files %s --%s " % (bamfilename, mapping)
    if mapping == "center":
        argstr += ("--nibble %s" % offset_or_nibble)
    else:
        argstr += ("--offset %s" % offset_or_nibble)

    mapdict = {
        "fiveprime" : FivePrimeMapFactory,
        "threeprime" : ThreePrimeMapFactory,
        "center" : CenterMapFactory,
        "fiveprime_variable" : VariableFivePrimeMapFactory,
    }
    
    args = alignment_file_parser.parse_args(shlex.split(argstr))
    ga = get_genome_array_from_args(args)
    #assert_equal(mapping,ga.map_fn.__mapping__)a

    assert_true(isinstance(ga.map_fn,mapdict[mapping]))

    # assure correct mapping function was applied
    # by comparing counts loaded from a bowtie file
    # under those map rules to a vector corresponding
    # to known answers for those rules
    assert_greater(ga.sum(),0)
    check_genome_array_against_vector(ga,mapping,offset_or_nibble)

def check_bowtie_mapping(mapping,offset_or_nibble):
    """Generator function to test import of bowtie files under many read mapping rules

    Parameters
    ----------
    mapping : str
        String specifying mapping rule. "fiveprime",
        "threeprime", "center", or "fiveprime_variable"

    offset_or_nibble : int or str
        Integer specifying offset or nibble (only for center mapping)
        or str specifying filename of variable offset dict text file
    """
    # assure correct mapping function was applied
    # by comparing counts loaded from a bowtie file
    # under those map rules to a vector corresponding
    # to known answers for those rules
    argstr = "--countfile_format bowtie --count_files %s --%s " % (bowtiefilename, mapping)
    if mapping == "center":
        argstr += ("--nibble %s" % offset_or_nibble)
    else:
        argstr += ("--offset %s" % offset_or_nibble)

    args = alignment_file_parser.parse_args(shlex.split(argstr))
    ga   = get_genome_array_from_args(args)
    assert_greater(ga.sum(),0)
    check_genome_array_against_vector(ga,mapping,offset_or_nibble)

@attr(test="unit")
def test_alignment_parser_bam_mapping():
    for mapping in ("fiveprime","threeprime","center"):
        yield check_bam_mapping, mapping, 0

    yield check_bam_mapping, "fiveprime_variable", offset_file

@attr(test="unit")
def test_alignment_parser_bowtie_mapping():
    for mapping in ("threeprime","fiveprime"):
        for offset in (0,15):
            yield check_bowtie_mapping, mapping, offset

    for nibble in (0,12):
        yield check_bowtie_mapping, "center", nibble 

    yield check_bowtie_mapping, "fiveprime_variable", offset_file

@attr(test="unit")
def test_alignment_parser_bam_fiveprime_variable_requires_offset():
    argstr = "--countfile_format BAM --count_files %s --fiveprime_variable" % bamfilename
    args = alignment_file_parser.parse_args(shlex.split(argstr))
    assert_raises(SystemExit,get_genome_array_from_args,args)

@attr(test="unit")
def test_alignment_parser_bowtie_fiveprime_varible_requires_offset():
    argstr = "--countfile_format bowtie --count_files %s --fiveprime_variable" % bowtiefilename
    args = alignment_file_parser.parse_args(shlex.split(argstr))
    assert_raises(SystemExit,get_genome_array_from_args,args)

def alignment_parser_check_normalize(fmt,filename):
    argstr = "--countfile_format %s --count_files %s --fiveprime --normalize" % (fmt,filename)
    args = alignment_file_parser.parse_args(shlex.split(argstr))
    ga = get_genome_array_from_args(args)
    assert_true(ga._normalize)

@attr(test="unit")
def test_alignment_parser_normalize():
    tups = [("wiggle",wigglefilename),
            ("BAM"   ,bamfilename),
            ("bowtie",bowtiefilename)]
    for fmt, filename in tups:
        yield alignment_parser_check_normalize, fmt, filename

@attr(test="unit")
def test_alignment_parser_input_choices():
    # kind of a hack. make sure input choices
    # are correct based upon whether or not they
    # appear in --help
    all_choices = ("bowtie","tagalign","wiggle","BAM")
    for chosen in itertools.chain.from_iterable(itertools.combinations(all_choices,K) for K in range(len(all_choices))):
        tmp_str = "  --countfile_format {%s}" % ",".join(list(chosen))
        parser = get_alignment_file_parser(input_choices=chosen)
        assert_in(tmp_str,parser.format_help())


#=============================================================================
# INDEX: tests for annotation file parsing
#=============================================================================

def check_annotation_input_choices(parser_fn):
    """Check that correct input choices are offered for various annotation file parsers
    and that GFF3-specific options are offered only when GFF3 is an option

    Parameters
    ----------
    parser_fn : function
        Function that returns an :py:class:`~argparse.ArgumentParser`
        of genome annotation files
    """
    all_choices = ["BED","BigBed","GTF2","GFF3","PSL"]
    for chosen in itertools.chain.from_iterable(itertools.combinations(all_choices,K) for K in range(len(all_choices))):
        tmp_str = "  --annotation_format {%s}" % (",".join(list(chosen)))
        parser = parser_fn(input_choices=chosen)
        assert_in(tmp_str,parser.format_help())
        
        for attr in annotation_file_gff_opts:
            args = parser.parse_args([])
            if "GFF3" in chosen:
                getattr(args,attr)
            else:
                assert_raises(AttributeError,getattr,args,attr)

@attr(test="unit")
def test_annotation_input_choices():
    pfuncs = [get_annotation_file_parser,
              get_segmentchain_file_parser,
             ]
    for parser_fn in pfuncs:
        yield check_annotation_input_choices, parser_fn

def check_get_ivcollection_from_args_return_type(parser,func,return_type,fmt):
    """Check return type of annotation file reading

    Parameters
    ----------
    parser : :py:class:`argparse.ArgumentParser`
        Argument parser for annotation files

    func : function
        Function that parsers the :py:class:`~argparse.ArgumentParser`
        into lists of accepted and rejected features of interest

    return_type : class
        Return type expected by function. |SegmentChain|, |Transcript|,
        or |BLAT_Alignment|

    fmt : str
        Format of input file ("BED", "BigBed", "GTF2", "GFF3", or "PSL")
    """
    argstr = "--annotation_format %s --annotation_files %s" % (fmt,MINI["%s_file" % fmt.lower()])
    args = parser.parse_args(shlex.split(argstr))
    transcripts = func(args)
    for tx in transcripts:
        assert_true(isinstance(tx,return_type))

@attr(test="unit")
def test_get_ivcollection_from_args_return_type():
    for fmt in ("BED","BigBed","GTF2","GFF3","PSL"):
        yield check_get_ivcollection_from_args_return_type, ivcollection_file_parser, get_segmentchains_from_args, SegmentChain, fmt

@attr(test="unit")
def test_get_transcript_from_args_return_type():
    for fmt in ("BED","BigBed","GTF2","GFF3"):
        yield check_get_ivcollection_from_args_return_type, annotation_file_parser, get_transcripts_from_args, Transcript, fmt

def check_get_transcript_from_args_add_three(parser,fmt,add_three,n=1):
    """Make sure stop codons are moved by three if ``add_three`` is specified
    in argument string and not moved by three if not specified.
    These are checked against expected values found in a BED file.
    
    Parameters
    ----------
    parser : :py:class:`argparse.ArgumentParser`
        Argument parser for annotation files

    fmt : str
        Format of input file ("BED", "BigBed", "GTF2", "GFF3", or "PSL")
        
    add_three : bool
        Whether or not to move stop codons by three nucleotides
    
    n : int
        Number of times input file was passed to parser (Default : 1)
    """
    # GTF file does not include stop codon; others do.
    # So, GTF file will be off by three relative to others
    # unless it includes explicit stop codon features
    #
    # we correct for that here
    gtf_offset = 3 if fmt == "GTF2" else 0
    bed_ref = MINI["bed_file"]
    dtmp = {}
    for line in open(bed_ref):
        items = line.strip("\n").split("\t")
        name = items[3]
        strand = items[5]
        cds_start = int(items[6])
        cds_end   = int(items[7])
        if cds_start < cds_end:
            if strand == "+":
                cds_end -= gtf_offset
            else:
                cds_start += gtf_offset
                
            if add_three == True:
                if strand == "+":
                    cds_end += 3
                else:
                    cds_start -= 3
            dtmp[name] = (cds_start,cds_end)
        else:
            dtmp[name] = (None,None)
    
    files = " ".join([MINI["%s_file" % fmt.lower()]] * n)
    argstr = "--annotation_format %s --annotation_files %s " % (fmt,files)
    if add_three == True:
        argstr += " --add_three"
    
    args = parser.parse_args(shlex.split(argstr))
    transcripts = list(get_transcripts_from_args(args))
    for tx in transcripts:
        name = tx.get_name()
        expected_start = dtmp[name][0]
        expected_end   = dtmp[name][1]
        found_start = tx.cds_genome_start
        found_end   = tx.cds_genome_end
        if found_start is None:
            assert_true(found_end is None)
            assert_true(expected_start is None,"Unequal CDS start (expected %s, got %s) on transcript %s" % (expected_start,found_start,name))
            assert_true(expected_end is None,"Unequal CDS stop (expected %s, got %s) on transcript %s" % (expected_end,found_end,name))
        else:
            assert_equal(found_start,expected_start,"Unequal CDS start (expected %s, got %s) on transcript %s" % (expected_start,found_start,name))
            assert_equal(found_end,expected_end,"Unequal CDS stop (expected %s, got %s) on transcript %s" % (expected_end,found_end,name))

    assert_equal(len(transcripts),n*len(dtmp),"Not all transcripts found in input. Expected %s. Got %s." % (n*len(dtmp),len(transcripts)))

@attr(test="unit")
def test_get_transcript_from_args_add_three():
    for fmt in ("BED","BigBed","GTF2","GFF3"):
        yield check_get_transcript_from_args_add_three, annotation_file_parser, fmt, False

    for fmt in ("BED","BigBed","GTF2","GFF3"):
        yield check_get_transcript_from_args_add_three, annotation_file_parser, fmt, True

@attr(test="unit")
def test_get_transcript_from_args_add_three_multiple_files():
    for fmt in ("BED","GTF2","GFF3"):
        yield check_get_transcript_from_args_add_three, annotation_file_parser, fmt, False, 3

@attr(test="unit")
def test_get_transcript_from_args_multiple_bigbed_raises_error():
    files = " ".join([MINI["bigbed_file"]] * 2)
    argstr = "--annotation_format BigBed --annotation_files %s " % (files)
    parser = get_annotation_file_parser()
    args = parser.parse_args(shlex.split(argstr))
    assert_raises(SystemExit,get_transcripts_from_args,args)

@attr(test="unit")
def check_get_transcript_from_args_add_three_tabix(fmt):
    # check endpoints of transcript and CDS
    # GTF2 and BED format
    bed_ref = REF_FILES["100transcripts_bed"]
    cds = {}
    tx  = {}
    for line in open(bed_ref):
        items = line.strip("\n").split("\t")
        name = items[3]

        tx_start = int(items[1])
        tx_end   = int(items[2])
        tx[name] = (tx_start,tx_end)
        
        cds_start = int(items[6])
        cds_end   = int(items[7])
        if cds_start < cds_end:
            cds[name] = (cds_start,cds_end)
        else:
            cds[name] = (None,None)
    
    ref_files = {
             "BED"  : REF_FILES["100transcripts_bed_tabix"],
             "GTF2" : REF_FILES["100transcripts_gtf_tabix"]
             }    
    parser = get_annotation_file_parser()
    argstr = "--annotation_format %s --annotation_files %s --tabix" % (fmt,ref_files[fmt])
    args = parser.parse_args(shlex.split(argstr))
    transcripts = list(get_transcripts_from_args(args))
    
    for my_tx in transcripts:
        name = my_tx.get_name()
        expected_cds_start, expected_cds_end = cds[name]
        
        found_cds_start = my_tx.cds_genome_start
        found_cds_end   = my_tx.cds_genome_end
        
        found_tx_start = my_tx.spanning_segment.start
        found_tx_end   = my_tx.spanning_segment.end
        
        expected_tx_start, expected_tx_end = tx[name]
        if found_cds_start is None:
            assert_true(found_cds_end is None)
            assert_true(expected_cds_start is None,"Unequal CDS start (expected %s, got %s) on transcript %s" % (expected_cds_start,found_cds_start,name))
            assert_true(expected_cds_end is None,"Unequal CDS stop (expected %s, got %s) on transcript %s" % (expected_cds_end,found_cds_end,name))
        else:
            assert_equal(found_cds_start,expected_cds_start,"Unequal CDS start (expected %s, got %s) on transcript %s" % (expected_cds_start,found_cds_start,name))
            assert_equal(found_cds_end,expected_cds_end,"Unequal CDS stop (expected %s, got %s) on transcript %s" % (expected_cds_end,found_cds_end,name))

        assert_equal(found_tx_start,expected_tx_start,"Unequal transcript start (expected %s, got %s) on transcript %s" % (expected_tx_start,found_tx_start,name))
        assert_equal(found_tx_end,expected_tx_end,"Unequal transcript end (expected %s, got %s) on transcript %s" % (expected_tx_end,found_tx_end,name))
    
        assert_equal(len(transcripts),len(cds),"Not all transcripts found in input. Expected %s. Got %s." % (len(cds),len(transcripts)))

@attr(test="unit")
def test_get_transcript_from_args_add_three_tabix():
    for fmt in ("BED","GTF2"):
        yield check_get_transcript_from_args_add_three_tabix, fmt

def check_arg_not_raises_error(name,func,args,params,kwargs):
    try:
        params = [args] + params
        func(*params,**kwargs)
        assert True
    except IOError: # because we're using dummy files
        pass

def check_arg_raises_error(name,func,args,params,kwargs):
    msg = "%s did not raise error with args '%s', params '%s' and kwargs '%s'" % (name,args,params,kwargs)
    params = [args] + params
    assert_raises(SystemExit,func,*params,**kwargs)

@attr(test="unit")
def test_annotation_not_sorted_raises_error_if_required():
    for func in get_transcripts_from_args, get_segmentchains_from_args:
        for fmt in ("BED","GTF2","GFF3"):
            argstr = "--annotation_format %s --annotation_files some_file" % fmt
            args = annotation_file_parser.parse_args(shlex.split(argstr))
            name = "require_sort_but_not_sorted_%s" % fmt
            yield check_arg_raises_error, name, func, args, [], {"require_sort" : True}

@attr(test="unit")
def test_annotation_sorted_raises_error_if_required():
    for func in get_transcripts_from_args, get_segmentchains_from_args:
        for fmt in ("BED","GTF2","GFF3"):
            argstr = "--sorted --annotation_format %s --annotation_files some_file" % fmt
            args = annotation_file_parser.parse_args(shlex.split(argstr))
            name = "require_sort_and_sorted_%s" % fmt
            yield check_arg_not_raises_error, name, func, args, [], {"require_sort" : True}
        for fmt in ("BED","GTF2","GFF3"):
            argstr = "--tabix --annotation_format %s --annotation_files some_file" % fmt
            args = annotation_file_parser.parse_args(shlex.split(argstr))
            name = "require_sort_and_tabix_%s" % fmt
            yield check_arg_not_raises_error, name, func, args, [], {"require_sort" : True}

@attr(test="unit")
def test_annotation_not_sorted_not_raises_error_if_not_required():
    for func in get_transcripts_from_args, get_segmentchains_from_args:
        for fmt in ("BED","GTF2","GFF3"):
            argstr = "--annotation_format %s --annotation_files some_file" % fmt
            args = annotation_file_parser.parse_args(shlex.split(argstr))
            name = "not_require_sort_and_not_sorted_%s" % fmt
            yield check_arg_not_raises_error, name, func, args, [], {"require_sort" : False}


#=============================================================================
# INDEX: tests for genome hash parsing
#=============================================================================

@attr(test="unit")
def test_mask_genome_hash_is_empty_if_no_file():
    argstr = "--mask_annotation_format BED"
    args = mask_file_parser.parse_args(shlex.split(argstr))
    mask_hash = get_genome_hash_from_mask_args(args)
    assert_true(isinstance(mask_hash,GenomeHash))
    assert_equal(len(mask_hash.feature_dict),0)

def check_mask_genome_hash(input_format,num_features):
    """Checks that a |GenomeHash| is correctly parsed from arguments,
    by making sure that it is of the right subclass, and that it
    contains the expected number of features

    Parameters
    ----------
    input_format : str
        Format of input file ("BED","BigBed", or "PSL")

    num_features : int
        Number of features in input file
    """
    argstr = "--mask_annotation_format %s --mask_annotation_files %s" % (input_format,MINI["%s_file" % input_format.lower()])
    args = mask_file_parser.parse_args(shlex.split(argstr))
    mask_hash = get_genome_hash_from_mask_args(args)

    # check correct return type and make sure features are populated
    if input_format == "BigBed":
        assert_true(isinstance(mask_hash,BigBedGenomeHash))
        assert_equal(mask_hash.bigbedreader.num_records,num_features)
    else:
        assert_true(isinstance(mask_hash,GenomeHash))
        assert_equal(len(mask_hash.feature_dict),num_features)

@attr(test="unit")
def test_mask_genome_hash():
    # psl is missing 2 features corresponding to introns in other files
    for input_format, num_features in [("BED",12),("BigBed",12),("PSL",10)]:
        yield check_mask_genome_hash, input_format, num_features



#=============================================================================
# INDEX: tests for sequence file parsing
#=============================================================================

def test_sequence_file_parser_opens_formats():
    for fmt,index in [("fasta",True),("fasta",False),("twobit",True)]:

        yield check_sequence_file_parser_opens_format, fmt, index

def check_sequence_file_parser_opens_format(fmt,index):
    argstr = "--sequence_file %s --sequence_format %s" % (REF_FILES["yeast_%s" % fmt],fmt)
    args = sequence_file_parser.parse_args(shlex.split(argstr))
    seqdict = get_seqdict_from_args(args,index)
    ref_chroms = sorted(TwoBitFile(REF_FILES["yeast_twobit"]).keys())
    msg = "Sequence file parser failed to open format '%s'." % fmt
    assert_list_equal(ref_chroms,sorted(seqdict.keys()),msg)


#=============================================================================
# INDEX: tests for PrefixNamespaceWrapper
#=============================================================================

@attr(test="unit")
class TestPrefixNamespaceWrapper(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # tuples of:
        #     argument name, argument type, default value
        cls.parser_options = [("foo",int,3),
                              ("bar",str,"barbarbar"),
                              ("baz",float,4.5),
                             ]

    def _get_parser(self,prefix=""):
        parser = argparse.ArgumentParser()
        for name,type_,val in self.parser_options:
            parser.add_argument("--%s%s" % (prefix,name),type=type_,default=val)

        return parser
    
    def _get_args(self,parser,fake_argv=[]):
        return parser.parse_args(fake_argv)

    def test_no_prefix(self):
        self._prefix_helper("")

    def test_prefix(self):
        self._prefix_helper("foofix_")

    def _prefix_helper(self,prefix):
        args = self._get_args(self._get_parser(prefix=prefix))
        wrapped = PrefixNamespaceWrapper(args,prefix)

        c = 0
        for name,type_,val in self.parser_options:
            c += 1
            wrapped_val = getattr(wrapped,name)
            raw_val = getattr(args,"%s%s" % (prefix,name))
            self.assertEqual(wrapped_val,val)
            self.assertEqual(raw_val,val)
            self.assertTrue(isinstance(wrapped_val,type_))

        self.assertEqual(c,len(self.parser_options))
