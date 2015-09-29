#!/usr/bin/env python
"""
"""
import numpy

from nose.plugins.attrib import attr
from nose.tools import assert_equal, assert_true, assert_set_equal, assert_less_equal
from numpy import nan

from plastid.util.io.filters import SkipBlankReader
from plastid.genomics.roitools import SegmentChain, GenomicSegment
from plastid.genomics.genome_hash import GenomeHash
from plastid.readers.gff import GFF3_TranscriptAssembler
from plastid.bin.metagene import window_landmark, \
                                     window_cds_start, \
                                     window_cds_stop, \
                                     group_regions_make_windows
from plastid.util.services.mini2to3 import cStringIO

_FLANKS = [(0,100),
           (50,100),
           (100,50),
           (100,0)
           ]
"""Regions upstream and downstream from landmarks to request"""

def check_equality(roi1,roi2,err_str=""):
    """Make sure two |SegmentChain| are equal by testing equality of position sets, chromosomes, and strands"""
    assert_set_equal(roi1.get_position_set(),roi2.get_position_set(),"%s ROIs unequal: unequal position sets %s vs %s" % (err_str,str(roi1),str(roi2)))
    if roi1.length > 0:
        assert_equal(roi1.spanning_segment.strand,roi2.spanning_segment.strand,"%s ROI strand mismatch: %s vs %s" % (err_str,roi1.spanning_segment.strand,roi2.spanning_segment.strand))
        assert_equal(roi1.spanning_segment.chrom,roi2.spanning_segment.chrom,"%s ROI chromosome mismatch: %s vs %s" % (err_str,roi1.spanning_segment.chrom,roi2.spanning_segment.chrom))

def check_window(tx_ivc,
                 known_roi,known_offset,known_ref_point,
                 flank_up,
                 flank_down,
                 test_method,
                 test_name,
                 ref_delta=0):
    """Helper function to test output of window landmark functions
    
    Parameters
    ----------
    tx_ivc : |SegmentChain|
        Test Transcript from which window will be derived
    
    known_roi : |SegmentChain|
        Reference output for ROI
    
    known_offset : int
        Known offset to start of ROI
    
    known_ref_point : (str,int,str) or numpy.nan
        Known offset to landmark in ROI as ("chromosome_name",position,"strand")
    
    flank_up : int
        Flank upstream of landmark to include in ROI
    
    flank_down : int
        Flank downstream of landmark to include in ROI
    
    test_method : function
        Function to test (e.g. :py:func:`window_cds_start`, py:func:`window_cds_stop`)
    
    test_name : str
        Name of test (for generating rich error output)
    
    ref_delta : int, optional
        Distance from reference landmark at which to center windows
    """
    err_str = ("Failed %s on %s (strand: '%s', up: %s, down: %s). " % (test_name,str(tx_ivc),tx_ivc.spanning_segment.strand,flank_up,flank_down)) + "%s unequal (%s vs %s)"
    test_roi, test_offset, test_ref_point = test_method(tx_ivc,flank_up,flank_down,ref_delta=ref_delta)
    check_equality(SegmentChain.from_str(known_roi),test_roi)

    # if no landmark
    if numpy.isnan(known_offset) or isinstance(known_ref_point,float) and numpy.isnan(known_ref_point):
        assert_true(numpy.isnan(test_offset),msg=err_str % ("offset",known_offset,test_offset))
        assert_true(numpy.isnan(test_ref_point),msg=err_str % ("ref_point",known_ref_point,test_ref_point))
    # if landmark
    else:
        assert_equal(known_offset,test_offset,msg=err_str % ("offset",known_offset,test_offset))
        assert_equal(known_ref_point,test_ref_point,msg=err_str % ("ref_point",known_ref_point,test_ref_point))

@attr(test="unit")
def test_window_cds_start():
    # test case, 1 transcript with CDS at first nucleotide
    # test case, 1 transcript with CDS internal
    # test case, 1 transcript with no CDS
    
    # for all, with and without changed upstream, downstream flanks
    for flank_up,flank_down in _FLANKS:
        for txid in _CDS_START_QUERIES:
            tx_ivc = _TRANSCRIPTS[txid]
            known_vals = _CDS_START_RESULTS["%s_%s_%s" % (tx_ivc.get_name(),flank_up,flank_down)]
            ltmp = [check_window,tx_ivc]
            ltmp.extend(known_vals)
            ltmp.append(flank_up)
            ltmp.append(flank_down)
            ltmp.append(window_cds_start)
            ltmp.append("window_cds_start")
            yield tuple(ltmp)

@attr(test="unit")
def test_window_cds_stop():
    # test case, 1 transcript with CDS stop at final nucleotide
    # test case, 1 transcript with CDS stop internal
    # test case, 1 transcript with no CDS
 
    # for all, with and without changed upstream, downstream flanks
    for flank_up,flank_down in _FLANKS:
        for txid in _CDS_STOP_QUERIES:
            tx_ivc = _TRANSCRIPTS[txid]
            known_vals = _CDS_STOP_RESULTS["%s_%s_%s" % (tx_ivc.get_name(),flank_up,flank_down)]
            ltmp = [check_window,tx_ivc]
            ltmp.extend(known_vals)
            ltmp.append(flank_up)
            ltmp.append(flank_down)
            ltmp.append(window_cds_stop)
            ltmp.append("window_cds_stop")
            yield tuple(ltmp)

@attr(test="unit")
def test_window_cds_stop_with_ref_delta():
    # test case, 1 transcript with CDS stop at final nucleotide
    # test case, 1 transcript with CDS stop internal
    # test case, 1 transcript with no CDS
 
    # for all, with and without changed upstream, downstream flanks
    for flank_up,flank_down in _FLANKS:
        for txid in _CDS_STOP_QUERIES:
            tx_ivc = _TRANSCRIPTS[txid]
            known_vals = _CDS_STOP_WITH_DELTA_RESULTS["%s_%s_%s" % (tx_ivc.get_name(),flank_up,flank_down)]
            ltmp = [check_window,tx_ivc]
            ltmp.extend(known_vals)
            ltmp.append(flank_up)
            ltmp.append(flank_down)
            ltmp.append(window_cds_stop)
            ltmp.append("window_cds_stop")
            ltmp.append(3)
            yield tuple(ltmp)

def check_window_landmark(ivc,landmark,flank_up=50,flank_down=100):
    test_roi, test_offset, (test_chrom,test_pos,test_strand) = window_landmark(ivc,flank_up,flank_down,landmark)
    
    # make sure position of interest matches what we gave window_landmark
    ref_chrom,ref_pos,ref_strand = ivc.get_genomic_coordinate(landmark)
    assert_equal(ref_chrom,test_chrom)
    assert_equal(ref_pos,test_pos)
    assert_equal(ref_strand,test_strand)

    # make sure position is in roi
    assert_true(test_pos in test_roi.get_position_set())
    
    # assure test roi length + offset == flank_up + flank_down
    # which it always should, unless the landmark is very close
    # to the edge of the transcript, in which case it should be smaller
    if landmark + flank_down <= test_roi.length:
        assert_equal(test_offset+test_roi.length,flank_up+flank_down)
    else:
        assert_less_equal(test_offset+test_roi.length,flank_up+flank_down)
    
    # test offset in roi is correct relative to offset and flank
    # only relevant for plus-strand
    roi_pos = test_roi.get_segmentchain_coordinate(test_chrom,test_pos,test_strand)
    assert_equal(roi_pos + test_offset,flank_up)

@attr(test="unit")
def test_window_landmark():
    # test cases: plus and minus-strand IVCs with splicing
    flank_up = 50
    flank_down = 100
    my_segmentchains = [SegmentChain(GenomicSegment("chrA",50,350,"+"),GenomicSegment("chrA",500,900,"+")),
               SegmentChain(GenomicSegment("chrA",50,350,"-"),GenomicSegment("chrA",500,900,"-")),
               ]
    for my_segmentchain in my_segmentchains:
        for landmark in range(0,700,50):
            yield check_window_landmark, my_segmentchain, landmark, flank_up, flank_down

def check_maximal_window(test_name,genome_hash,test_group,result_groups,flank_up,flank_down):
    """
    test_name : str
        Descriptive name of test

    genome_hash : GenomeHash
        Mask hash

    test_group : list
        List of transcript IDs, referring to transcripts in the GFF text above

    result_groups : list
        list of tuples of (region_str, aligment_offset, window_length) expected from 
        maximal spanning window output

    flank_up : int
        Bases to include upstream of landmark in window

    flank_down : int
        bases to include downstream of landmark in window
    """
    # table keys:
    #    gene_id
    #    window_size
    #    roi
    #    masked
    #    alignment_offset
    #    zero_point
    err_str = ("Failed %s (up: %s, down: %s). " % (test_name,flank_up,flank_down)) + "%s unequal (%s vs %s)"
    tx_ivcs = (_TRANSCRIPTS[X] for X in test_group)
    roi_table = group_regions_make_windows(tx_ivcs,genome_hash,flank_up,flank_down,window_cds_start)
    roi_table.sort(columns=["region"],inplace=True)
    trows = [X[1] for X in roi_table.iterrows()]
    result_groups = sorted(result_groups,key=lambda x: x[0])
    REGION = 0

    c = 0
    
    for n,result_group in enumerate(result_groups):
        # if no landmark
        if numpy.isnan(result_group[1]) or numpy.isnan(result_group[2]):
            c += 1 # increment counter for input that will have no output
            
        # if landmark
        else:
            check_equality(SegmentChain.from_str(result_group[0]),
                           SegmentChain.from_str(trows[n-c]["region"]),test_name) 
            assert_equal(result_group[1],
                         trows[n-c]["alignment_offset"],
                         msg=err_str % ("offset",result_group[1],trows[n-c]["alignment_offset"]))
            assert_equal(result_group[2],
                         trows[n-c]["zero_point"],
                         msg=err_str % ("ref_point",result_group[1],trows[n-c]["zero_point"]))
            if len(result_group) == 4:
                assert_equal(result_group[3],
                             trows[n-c]["masked"],
                             msg=err_str % ("mask",result_group[3],trows[n-c]["masked"]))

    assert_equal(n+1-c,len(roi_table))
    
@attr(test="unit")
def test_group_regions_make_windows_finds_maximal_window():
    # test case, 3 transcripts same start
    # test case, 3 transcripts, different starts
    # test case, 3 transcripts, 2 share start, 1 different
    # test case including ncRNAs
    # for all, with and without changed upstream, downstream flanks
    empty_hash = GenomeHash([])
    for flank_up,flank_down in _FLANKS:
        for test_name, test_group in _DO_GENERATE_MAX_WINDOW.items():
            result_group = _DO_GENERATE_MAX_WINDOW_RESULTS["%s_%s_%s" % (test_name,flank_up,flank_down)]
            yield check_maximal_window, test_name, empty_hash, test_group, [result_group], flank_up, flank_down

@attr(test="unit")
def test_group_regions_make_windows_multiple_genes():
    empty_hash = GenomeHash([])
    flank_up = 50
    flank_down = 100
    for test_name, test_group in _DO_GENERATE_MULTI_GENE.items():
        result_groups = _DO_GENERATE_MULTI_GENE_RESULTS["%s_%s_%s" % (test_name,flank_up,flank_down)]
        yield check_maximal_window, test_name, empty_hash, test_group, result_groups, flank_up, flank_down

@attr(test="unit")
def test_group_regions_make_windows_finds_masks():
    crossmap = GenomeHash(_MASKS)
    for flank_up,flank_down in _FLANKS:
        for test_name, test_group in _DO_GENERATE_MAX_WINDOW.items():
            result_group = _DO_GENERATE_MAX_WINDOW_RESULTS_MASKED["%s_%s_%s" % (test_name,flank_up,flank_down)]
            yield check_maximal_window, test_name, crossmap, test_group, [result_group], flank_up, flank_down


#===============================================================================
# INDEX: test data
#===============================================================================

_MASKS = [
    SegmentChain.from_str("2L:7985694-7985744(+)"),         
    SegmentChain.from_str("3R:4519879-4519891(-)"),
    SegmentChain.from_str("4:50-50000(+)"),
]


_TRANSCRIPTS_GFF = """##gff-version 3
3R    FlyBase    mRNA    4517211    4523544    .    -    .    ID=FBtr0081950;Name=hb-RB;Parent=FBgn0001180;Alias=FBtr0002097,FBtr0002098,CG9786-RB,hb[+]R2.8;Dbxref=FlyBase_Annotation_IDs:CG9786-RB,REFSEQ:NM_169234;score_text=Strongly Supported;score=11
3R    FlyBase    exon    4517211    4519894    .    -    .    Name=hb:2;Parent=FBtr0081950;parent_type=mRNA
3R    FlyBase    CDS    4517600    4519876    .    -    0    Name=hb-cds;Parent=FBtr0081950;parent_type=mRNA
3R    FlyBase    exon    4523048    4523544    .    -    .    Name=hb:4;Parent=FBtr0081950;parent_type=mRNA

3R    FlyBase    mRNA    4516702    4520322    .    -    .    ID=FBtr0081951;Name=hb-RA;Parent=FBgn0001180;Alias=FBtr0002096,FBtr0002097,CG9786-RA,hb[+]R3.2;Dbxref=FlyBase_Annotation_IDs:CG9786-RA,REFSEQ:NM_169233;score_text=Strongly Supported;score=11
3R    FlyBase    exon    4516702    4519894    .    -    .    Name=hb:1;Parent=FBtr0081951;parent_type=mRNA
3R    FlyBase    CDS    4517600    4519876    .    -    0    Name=hb-cds;Parent=FBtr0081951;parent_type=mRNA
3R    FlyBase    exon    4520178    4520322    .    -    .    Name=hb:3;Parent=FBtr0081951;parent_type=mRNA

3R    FlyBase    mRNA    4516702    4520322    .    -    .    ID=FBtr0081951_alt_start;Name=hb-RA_alt_start;Parent=FBgn0001180;Alias=FBtr0002096,FBtr0002097,CG9786-RA,hb[+]R3.2;Dbxref=FlyBase_Annotation_IDs:CG9786-RA,REFSEQ:NM_169233;score_text=Strongly Supported;score=11
3R    FlyBase    exon    4516702    4519894    .    -    .    Name=hb:1;Parent=FBtr0081951_alt_start;parent_type=mRNA
3R    FlyBase    CDS     4517600    4519816    .    -    0    Name=hb-cds;Parent=FBtr0081951_alt_start;parent_type=mRNA
3R    FlyBase    exon    4520178    4520322    .    -    .    Name=hb:3;Parent=FBtr0081951_alt_start;parent_type=mRNA

3R    FlyBase    mRNA    4517211    4523544    .    -    .    ID=FBtr0081950_no_cds;Name=hb-RB_no_cds;Parent=FBgn0001180;Alias=FBtr0002097,FBtr0002098,CG9786-RB,hb[+]R2.8;Dbxref=FlyBase_Annotation_IDs:CG9786-RB,REFSEQ:NM_169234;score_text=Strongly Supported;score=11
3R    FlyBase    exon    4517211    4519894    .    -    .    Name=hb_no_cds:2;Parent=FBtr0081950_no_cds;parent_type=mRNA
3R    FlyBase    exon    4523048    4523544    .    -    .    Name=hb_no_cds:4;Parent=FBtr0081950_no_cds;parent_type=mRNA

3R    FlyBase    mRNA    4517590    4523544    .    -    .    ID=FBtr0081950_short_utr;Name=hb-RB_short_utr;Parent=FBgn0001180;Alias=FBtr0002097,FBtr0002098,CG9786-RB,hb[+]R2.8;Dbxref=FlyBase_Annotation_IDs:CG9786-RB,REFSEQ:NM_169234;score_text=Strongly Supported;score=11
3R    FlyBase    exon    4517590    4519894    .    -    .    Name=hb:2;Parent=FBtr0081950_short_utr;parent_type=mRNA
3R    FlyBase    CDS    4517600    4519876    .    -    0    Name=hb-cds;Parent=FBtr0081950_short_utr;parent_type=mRNA

3R    FlyBase    mRNA    4517600    4519876    .    -    .    ID=FBtr0081950_no_utr;Name=hb-RB_no_utr;Parent=FBgn0001180;Alias=FBtr0002097,FBtr0002098,CG9786-RB,hb[+]R2.8;Dbxref=FlyBase_Annotation_IDs:CG9786-RB,REFSEQ:NM_169234;score_text=Strongly Supported;score=11
3R    FlyBase    exon    4517600    4519876    .    -    .    Name=hb:2;Parent=FBtr0081950_no_utr;parent_type=mRNA
3R    FlyBase    CDS    4517600    4519876    .    -    0    Name=hb-cds;Parent=FBtr0081950_no_utr;parent_type=mRNA

3R    FlyBase    mRNA    4517211    4523544    .    -    .    ID=FBtr0081950_at_splice;Name=hb-RB_at_splice;Parent=FBgn0001180;Alias=FBtr0002097,FBtr0002098,CG9786-RB,hb[+]R2.8;Dbxref=FlyBase_Annotation_IDs:CG9786-RB,REFSEQ:NM_169234;score_text=Strongly Supported;score=11
3R    FlyBase    exon    4523500    4523544    .    -    .    Name=hb:3;Parent=FBtr0081950_at_splice;parent_type=mRNA
3R    FlyBase    exon    4517211    4517241    .    -    .    Name=hb:2;Parent=FBtr0081950_at_splice;parent_type=mRNA
3R    FlyBase    exon    4517600    4519876    .    -    .    Name=hb:2g;Parent=FBtr0081950_at_splice;parent_type=mRNA
3R    FlyBase    CDS    4517600    4519876    .    -    0    Name=hb-cds;Parent=FBtr0081950_at_splice;parent_type=mRNA

2L    FlyBase    mRNA    7983946    7986785    .    +    .    ID=FBtr0300185;Name=CG7231-RE;Parent=FBgn0031968;Dbxref=FlyBase_Annotation_IDs:CG7231-RE,REFSEQ:NM_164776;score_text=Moderately Supported;score=7
2L    FlyBase    CDS    7984768    7984810    .    +    0    Name=CG7231-cds;Parent=FBtr0300185;parent_type=mRNA
2L    FlyBase    CDS    7985553    7985768    .    +    2    Name=CG7231-cds;Parent=FBtr0300185;parent_type=mRNA
2L    FlyBase    CDS    7986139    7986705    .    +    0    Name-=CG7231-cds;Parent=FBtr0300185;parent_type=mRNA
2L    FlyBase    CDS    7985834    7986084    .    +    2    Name=CG7231-cds;Parent=FBtr0300185;parent_type=mRNA
2L    FlyBase    exon    7983946    7984238    .    +    .    Name=CG7231:1;Parent=FBtr0300185;parent_type=mRNA
2L    FlyBase    exon    7984541    7984810    .    +    .    Name=CG7231:2;Parent=FBtr0300185;parent_type=mRNA
2L    FlyBase    exon    7985553    7985768    .    +    .    Name=CG7231:6b;Parent=FBtr0300185;parent_type=mRNA
2L    FlyBase    exon    7985834    7986084    .    +    .    Name=CG7231:8;Parent=FBtr0300185;parent_type=mRNA
2L    FlyBase    exon    7986139    7986785    .    +    .    Name=CG7231:9b;Parent=FBtr0300185;parent_type=mRNA

2L    FlyBase    mRNA    7985283    7986785    .    +    .    ID=FBtr0079531;Name=CG7231-RC;Parent=FBgn0031968;Dbxref=FlyBase_Annotation_IDs:CG7231-RC,REFSEQ:NM_205923;score_text=Moderately Supported;score=7
2L    FlyBase    exon    7985283    7985433    .    +    .    Name=CG7231:4;Parent=FBtr0079531;parent_type=mRNA
2L    FlyBase    CDS    7985675    7985768    .    +    0    Name=CG7231-cdsc;Parent=FBtr0079531;parent_type=mRNA
2L    FlyBase    CDS    7985834    7986075    .    +    2    Name=CG7231-cdsc;Parent=FBtr0079531;parent_type=mRNA
2L    FlyBase    CDS    7986139    7986705    .    +    0    Name=CG7231-cdsc;Parent=FBtr0079531;parent_type=mRNA
2L    FlyBase    exon    7985553    7985768    .    +    .    Name=CG7231:6c;Parent=FBtr0079531;parent_type=mRNA
2L    FlyBase    exon    7985834    7986075    .    +    .    Name=CG7231:7c;Parent=FBtr0079531;parent_type=mRNA
2L    FlyBase    exon    7986139    7986785    .    +    .    Name=CG7231:9c;Parent=FBtr0079531;parent_type=mRNA

2L    FlyBase    mRNA    7985283    7986785    .    +    .    ID=FBtr0306336;Name=CG7231-RF;Parent=FBgn0031968;Dbxref=REFSEQ:NM_001259004,FlyBase_Annotation_IDs:CG7231-RF;score_text=Moderately Supported;score=7
2L    FlyBase    exon    7985283    7985441    .    +    .    Name=CG7231:5;Parent=FBtr0306336;parent_type=mRNA
2L    FlyBase    CDS    7985675    7985768    .    +    0    Name=CG7231-cds;Parent=FBtr0306336;parent_type=mRNA
2L    FlyBase    CDS    7985834    7986075    .    +    2    Name=CG7231-cds;Parent=FBtr0306336;parent_type=mRNA
2L    FlyBase    CDS    7986139    7986705    .    +    0    Name=CG7231-cds;Parent=FBtr0306336;parent_type=mRNA
2L    FlyBase    exon    7985553    7985768    .    +    .    Name=CG7231:6;Parent=FBtr0306336;parent_type=mRNA
2L    FlyBase    exon    7985834    7986075    .    +    .    Name=CG7231:7;Parent=FBtr0306336;parent_type=mRNA
2L    FlyBase    exon    7986139    7986785    .    +    .    Name=CG7231:9;Parent=FBtr0306336;parent_type=mRNA

2L    FlyBase    mRNA    7985665    7986715    .    +    .    ID=FBtr0079531_short_utr;Name=CG7231-RC;Parent=FBgn0031968;Dbxref=FlyBase_Annotation_IDs:CG7231-RC,REFSEQ:NM_205923;score_text=Moderately Supported;score=7
2L    FlyBase    exon    7985665    7985768    .    +    .    Name=CG7231:6;Parent=FBtr0079531_short_utr;parent_type=mRNA
2L    FlyBase    exon    7985834    7986075    .    +    .    Name=CG7231:7;Parent=FBtr0079531_short_utr;parent_type=mRNA
2L    FlyBase    exon    7986139    7986715    .    +    .    Name=CG7231:9;Parent=FBtr0079531_short_utr;parent_type=mRNA
2L    FlyBase    CDS    7985675    7985768    .    +    0    Name=CG7231-cds;Parent=FBtr0079531_short_utr;parent_type=mRNA
2L    FlyBase    CDS    7985834    7986075    .    +    2    Name=CG7231-cds;Parent=FBtr0079531_short_utr;parent_type=mRNA
2L    FlyBase    CDS    7986139    7986705    .    +    0    Name=CG7231-cds;Parent=FBtr0079531_short_utr;parent_type=mRNA

2L    FlyBase    mRNA    7985283    7986785    .    +    .    ID=FBtr0079531_no_cds;Name=CG7231-RC;Parent=FBgn0031968;Dbxref=FlyBase_Annotation_IDs:CG7231-RC,REFSEQ:NM_205923;score_text=Moderately Supported;score=7
2L    FlyBase    exon    7985283    7985433    .    +    .    Name=CG7231:4;Parent=FBtr0079531_no_cds;parent_type=mRNA
2L    FlyBase    exon    7985553    7985768    .    +    .    Name=CG7231:6;Parent=FBtr0079531_no_cds;parent_type=mRNA
2L    FlyBase    exon    7985834    7986075    .    +    .    Name=CG7231:7;Parent=FBtr0079531_no_cds;parent_type=mRNA
2L    FlyBase    exon    7986139    7986785    .    +    .    Name=CG7231:9;Parent=FBtr0079531_no_cds;parent_type=mRNA

2L    FlyBase    mRNA    7985675    7986075    .    +    .    ID=FBtr0079531_no_utr;Name=CG7231-RC;Parent=FBgn0031968;Dbxref=FlyBase_Annotation_IDs:CG7231-RC,REFSEQ:NM_205923;score_text=Moderately Supported;score=7
2L    FlyBase    exon    7985675    7985768    .    +    .    Name=CG7231:6;Parent=FBtr0079531_no_utr;parent_type=mRNA
2L    FlyBase    exon    7985834    7986075    .    +    .    Name=CG7231:7;Parent=FBtr0079531_no_utr;parent_type=mRNA
2L    FlyBase    CDS    7985675    7985768    .    +    0    Name=CG7231-cds;Parent=FBtr0079531_no_utr;parent_type=mRNA
2L    FlyBase    CDS    7985834    7986075    .    +    2    Name=CG7231-cds;Parent=FBtr0079531_no_utr;parent_type=mRNA
2L    FlyBase    CDS    7986139    7986705    .    +    0    Name=CG7231-cds;Parent=FBtr0079531_no_utr;parent_type=mRNA

2L    FlyBase    mRNA    7985283    7986785    .    +    .    ID=FBtr0079531_at_splice;Name=CG7231-RC;Parent=FBgn0031968;Dbxref=FlyBase_Annotation_IDs:CG7231-RC,REFSEQ:NM_205923;score_text=Moderately Supported;score=7
2L    FlyBase    exon    7985283    7985433    .    +    .    Name=CG7231:4;Parent=FBtr0079531_at_splice;parent_type=mRNA
2L    FlyBase    exon    7985675    7985768    .    +    .    Name=CG7231:6;Parent=FBtr0079531_at_splice;parent_type=mRNA
2L    FlyBase    exon    7985834    7986075    .    +    .    Name=CG7231:7;Parent=FBtr0079531_at_splice;parent_type=mRNA
2L    FlyBase    exon    7986139    7986705    .    +    .    Name=CG7231:9;Parent=FBtr0079531_at_splice;parent_type=mRNA
2L    FlyBase    exon    7986775    7986785    .    +    .    Name=CG7231:9;Parent=FBtr0079531_at_splice;parent_type=mRNA
2L    FlyBase    CDS    7985675    7985768    .    +    0    Name=CG7231-cds;Parent=FBtr0079531_at_splice;parent_type=mRNA
2L    FlyBase    CDS    7985834    7986075    .    +    2    Name=CG7231-cds;Parent=FBtr0079531_at_splice;parent_type=mRNA
2L    FlyBase    CDS    7986139    7986705    .    +    0    Name=CG7231-cds;Parent=FBtr0079531_at_splice;parent_type=mRNA

2L    FlyBase    mRNA    8997641    8998544    .    +    .    ID=FBtr0079763;Name=CG18661-RA;Parent=FBgn0040964;Alias=CG18661-RB;Dbxref=FlyBase_Annotation_IDs:CG18661-RA,REFSEQ:NM_144305;score_text=Strongly Supported;score=15
2L    FlyBase    mRNA    8997641    8998544    .    +    .    ID=FBtr0303880;Name=CG18661-RB;Parent=FBgn0040964;Dbxref=FlyBase_Annotation_IDs:CG18661-RB,REFSEQ:NM_001201814;score_text=Strongly Supported;score=15
2L    FlyBase    exon    8997641    8997668    .    +    .    ID=FBgn0040964:1;Name=CG18661:1;Parent=FBtr0079763,FBtr0303880;parent_type=mRNA
2L    FlyBase    CDS    8997664    8997668    .    +    0    ID=CDS_FBgn0040964:1_866;Name=CG18661-cds;Parent=FBtr0303880;parent_type=mRNA
2L    FlyBase    CDS    8997664    8997668    .    +    0    ID=CDS_FBgn0040964:1_867;Name=CG18661-cds;Parent=FBtr0079763;parent_type=mRNA
2L    FlyBase    exon    8997723    8998544    .    +    .    ID=FBgn0040964:2;Name=CG18661:2;Parent=FBtr0079763;parent_type=mRNA
2L    FlyBase    CDS    8997723    8998386    .    +    1    ID=CDS_FBgn0040964:2_867;Name=CG18661-cds;Parent=FBtr0079763;parent_type=mRNA
2L    FlyBase    exon    8997762    8998544    .    +    .    ID=FBgn0040964:3;Name=CG18661:3;Parent=FBtr0303880;parent_type=mRNA
2L    FlyBase    CDS    8997762    8998386    .    +    1    ID=CDS_FBgn0040964:3_866;Name=CG18661-cds;Parent=FBtr0303880;parent_type=mRNA

2L    FlyBase    mRNA    9397110    9397988    .    -    .    ID=FBtr0079813;Name=CG4438-RA;Parent=FBgn0032115;
2L    FlyBase    mRNA    9397110    9397988    .    -    .    ID=FBtr0303900;Name=CG4438-RB;Parent=FBgn0032115;
2L    FlyBase    exon    9397110    9397772    .    -    .    ID=FBgn0032115:2;Name=CG4438:2;Parent=FBtr0079813;parent_type=mRNA
2L    FlyBase    exon    9397110    9397764    .    -    .    ID=FBgn0032115:1;Name=CG4438:1;Parent=FBtr0303900;parent_type=mRNA
2L    FlyBase    CDS    9397177    9397746    .    -    0    ID=CDS_FBgn0032115:1_867;Name=CG4438-cds;Parent=FBtr0303900;parent_type=mRNA
2L    FlyBase    CDS    9397177    9397746    .    -    0    ID=CDS_FBgn0032115:2_867;Name=CG4438-cds;Parent=FBtr0079813;parent_type=mRNA
2L    FlyBase    exon    9397833    9397988    .    -    .    ID=FBgn0032115:3;Name=CG4438:3;Parent=FBtr0079813,FBtr0303900;parent_type=mRNA
""".replace("    ","\t")
"""GFF of transcripts used in these tests"""

_TRANSCRIPTS = { X.get_name() : X for X in GFF3_TranscriptAssembler(SkipBlankReader(cStringIO.StringIO(_TRANSCRIPTS_GFF))) }
"""|Transcript| representation of transcripts used in these tests"""

_CDS_START_QUERIES = [ 
    "FBtr0081950",           # minus-strand, CDS internal, UTR longer than flank
    "FBtr0081950_short_utr", # minus-strand, UTR will be shorter than flank
    "FBtr0081950_no_cds",    # minus-strand, no CDS
    "FBtr0081950_no_utr",    # minus-strand, CDS start at end of tx
    "FBtr0081950_at_splice", # minus-strand, CDS start at splice junction   
    "FBtr0079531",                  
    "FBtr0079531_short_utr",                  
    "FBtr0079531_no_cds",                  
    "FBtr0079531_no_utr",                  
    "FBtr0079531_at_splice",    
]
"""IDs of individual |Transcript| s to use in :py:func:`test_window_cds_start`"""

_CDS_START_RESULTS = {
    'FBtr0081950_0_100'           : ('3R:4519776-4519876(-)', 0, ('3R', 4519875, '-')),
    'FBtr0081950_100_0'           : ('3R:4519876-4519894^4523047-4523129(-)',0,('3R', 4519875, '-')),
    'FBtr0081950_100_50'          : ('3R:4519826-4519894^4523047-4523129(-)',0,('3R', 4519875, '-')),
    'FBtr0081950_50_100'          : ('3R:4519776-4519894^4523047-4523079(-)',0,('3R', 4519875, '-')),

    'FBtr0081950_no_cds_0_100'    : ('na', nan, nan),
    'FBtr0081950_no_cds_100_0'    : ('na', nan, nan),
    'FBtr0081950_no_cds_100_50'   : ('na', nan, nan),
    'FBtr0081950_no_cds_50_100'   : ('na', nan, nan),
 
    'FBtr0081950_no_utr_0_100'    : ('3R:4519776-4519876(-)',0,('3R', 4519875, '-')),
    'FBtr0081950_no_utr_100_0'    : ('na', 100, ('3R', 4519875, '-')),
    'FBtr0081950_no_utr_100_50'   : ('3R:4519826-4519876(-)',100,('3R', 4519875, '-')),
    'FBtr0081950_no_utr_50_100'   : ('3R:4519776-4519876(-)',50,('3R', 4519875, '-')),

    'FBtr0081950_short_utr_0_100' : ('3R:4519776-4519876(-)',0,('3R', 4519875, '-')),
    'FBtr0081950_short_utr_100_0' : ('3R:4519876-4519894(-)',82,('3R', 4519875, '-')),
    'FBtr0081950_short_utr_100_50': ('3R:4519826-4519894(-)',82,('3R', 4519875, '-')),
    'FBtr0081950_short_utr_50_100': ('3R:4519776-4519894(-)',32,('3R', 4519875, '-')),

    'FBtr0081950_at_splice_0_100': ('3R:4519776-4519876(-)',0,('3R', 4519875, '-')),
    'FBtr0081950_at_splice_100_0': ('3R:4523499-4523544(-)',55,('3R', 4519875, '-')),
    'FBtr0081950_at_splice_100_50': ('3R:4519826-4519876^4523499-4523544(-)',55,('3R', 4519875, '-')),
    'FBtr0081950_at_splice_50_100': ('3R:4519776-4519876^4523499-4523544(-)',5,('3R', 4519875, '-')),    

    'FBtr0081950_no_utr_100_0': ('na', 100, ('3R', 4519875, '-')),
    'FBtr0081950_no_utr_100_50': ('3R:4519826-4519876(-)',100,('3R', 4519875, '-')),
    'FBtr0081950_no_utr_50_100': ('3R:4519776-4519876(-)',50,('3R', 4519875, '-')),

    'FBtr0081950_short_utr_0_100': ('3R:4519776-4519876(-)',0,('3R', 4519875, '-')),
    'FBtr0081950_short_utr_100_0': ('3R:4519876-4519894(-)',82,('3R', 4519875, '-')),
    'FBtr0081950_short_utr_100_50': ('3R:4519826-4519894(-)',82,('3R', 4519875, '-')),
    'FBtr0081950_short_utr_50_100': ('3R:4519776-4519894(-)',32,('3R', 4519875, '-')),
    
    'FBtr0079531_0_100': ('2L:7985674-7985768^7985833-7985839(+)',0,('2L', 7985674, '+')),
    'FBtr0079531_100_0': ('2L:7985574-7985674(+)', 0, ('2L', 7985674, '+')),
    'FBtr0079531_100_50': ('2L:7985574-7985724(+)', 0, ('2L', 7985674, '+')),
    'FBtr0079531_50_100': ('2L:7985624-7985768^7985833-7985839(+)',0,('2L', 7985674, '+')),
    'FBtr0079531_at_splice_0_100': ('2L:7985674-7985768^7985833-7985839(+)',0,('2L', 7985674, '+')),
    'FBtr0079531_at_splice_100_0': ('2L:7985333-7985433(+)',0,('2L', 7985674, '+')),
    'FBtr0079531_at_splice_100_50': ('2L:7985333-7985433^7985674-7985724(+)',0,('2L', 7985674, '+')),
    'FBtr0079531_at_splice_50_100': ('2L:7985383-7985433^7985674-7985768^7985833-7985839(+)',0,('2L', 7985674, '+')),
    'FBtr0079531_no_cds_0_100': ('na', nan, nan),
    'FBtr0079531_no_cds_100_0': ('na', nan, nan),
    'FBtr0079531_no_cds_100_50': ('na', nan, nan),
    'FBtr0079531_no_cds_50_100': ('na', nan, nan),
    'FBtr0079531_no_utr_0_100': ('2L:7985674-7985768^7985833-7985839(+)',0,('2L', 7985674, '+')),
    'FBtr0079531_no_utr_100_0': ('na', 100, ('2L', 7985674, '+')),
    'FBtr0079531_no_utr_100_50': ('2L:7985674-7985724(+)',100,('2L', 7985674, '+')),
    'FBtr0079531_no_utr_50_100': ('2L:7985674-7985768^7985833-7985839(+)',50,('2L', 7985674, '+')),
    'FBtr0079531_short_utr_0_100': ('2L:7985674-7985768^7985833-7985839(+)',0,('2L', 7985674, '+')),
    'FBtr0079531_short_utr_100_0': ('2L:7985664-7985674(+)',90,('2L', 7985674, '+')),
    'FBtr0079531_short_utr_100_50': ('2L:7985664-7985724(+)',90,('2L', 7985674, '+')),
    'FBtr0079531_short_utr_50_100': ('2L:7985664-7985768^7985833-7985839(+)',40,('2L', 7985674, '+')),     
}
"""Expected results of CDS start queries"""

_CDS_STOP_QUERIES = [ 
    "FBtr0081950",           # minus-strand, CDS internal, UTR longer than flank
    "FBtr0081950_short_utr", # minus-strand, UTR will be shorter than flank
    "FBtr0081950_no_cds",    # minus-strand, no CDS
    "FBtr0081950_no_utr",    # minus-strand, CDS start at end of tx
    "FBtr0081950_at_splice", # minus-strand, CDS start at splice junction   
    "FBtr0079531",                  
    "FBtr0079531_short_utr",                  
    "FBtr0079531_no_cds",                  
    "FBtr0079531_no_utr",                  
    "FBtr0079531_at_splice",    
]
"""IDs of individual |Transcript| s to use in :py:func:`test_window_cds_stop`"""

_CDS_STOP_RESULTS = {
    'FBtr0079531_0_100': ('2L:7986702-7986785(+)', 0, ('2L', 7986702, '+')),
    'FBtr0079531_100_0': ('2L:7986602-7986702(+)', 0, ('2L', 7986702, '+')),
    'FBtr0079531_100_50': ('2L:7986602-7986752(+)', 0, ('2L', 7986702, '+')),
    'FBtr0079531_50_100': ('2L:7986652-7986785(+)', 0, ('2L', 7986702, '+')),
    'FBtr0079531_at_splice_0_100': ('2L:7986702-7986705^7986774-7986785(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_at_splice_100_0': ('2L:7986602-7986702(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_at_splice_100_50': ('2L:7986602-7986705^7986774-7986785(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_at_splice_50_100': ('2L:7986652-7986705^7986774-7986785(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_no_cds_0_100': ('na', nan, nan),
    'FBtr0079531_no_cds_100_0': ('na', nan, nan),
    'FBtr0079531_no_cds_100_50': ('na', nan, nan),
    'FBtr0079531_no_cds_50_100': ('na', nan, nan),
    'FBtr0079531_no_utr_0_100': ('2L:7986702-7986705(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_no_utr_100_0': ('2L:7986602-7986702(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_no_utr_100_50': ('2L:7986602-7986705(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_no_utr_50_100': ('2L:7986652-7986705(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_short_utr_0_100': ('2L:7986702-7986715(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_short_utr_100_0': ('2L:7986602-7986702(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_short_utr_100_50': ('2L:7986602-7986715(+)',0,('2L', 7986702, '+')),
    'FBtr0079531_short_utr_50_100': ('2L:7986652-7986715(+)',0,('2L', 7986702, '+')),
    'FBtr0081950_0_100': ('3R:4517502-4517602(-)', 0, ('3R', 4517601, '-')),
    'FBtr0081950_100_0': ('3R:4517602-4517702(-)', 0, ('3R', 4517601, '-')),
    'FBtr0081950_100_50': ('3R:4517552-4517702(-)', 0, ('3R', 4517601, '-')),
    'FBtr0081950_50_100': ('3R:4517502-4517652(-)', 0, ('3R', 4517601, '-')),
    'FBtr0081950_at_splice_0_100': ('3R:4517210-4517241^4517599-4517602(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_at_splice_100_0': ('3R:4517602-4517702(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_at_splice_100_50': ('3R:4517210-4517241^4517599-4517702(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_at_splice_50_100': ('3R:4517210-4517241^4517599-4517652(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_no_cds_0_100': ('na', nan, nan),
    'FBtr0081950_no_cds_100_0': ('na', nan, nan),
    'FBtr0081950_no_cds_100_50': ('na', nan, nan),
    'FBtr0081950_no_cds_50_100': ('na', nan, nan),
    'FBtr0081950_no_utr_0_100': ('3R:4517599-4517602(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_no_utr_100_0': ('3R:4517602-4517702(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_no_utr_100_50': ('3R:4517599-4517702(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_no_utr_50_100': ('3R:4517599-4517652(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_short_utr_0_100': ('3R:4517589-4517602(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_short_utr_100_0': ('3R:4517602-4517702(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_short_utr_100_50': ('3R:4517589-4517702(-)',0,('3R', 4517601, '-')),
    'FBtr0081950_short_utr_50_100': ('3R:4517589-4517652(-)',0,('3R', 4517601, '-'))             
}
"""Expected results of queries"""

_CDS_STOP_WITH_DELTA_RESULTS = {
    'FBtr0079531_0_100': ('2L:7986705-7986785(+)', 0, ('2L', 7986705, '+')),
    'FBtr0079531_100_0': ('2L:7986605-7986705(+)', 0, ('2L', 7986705, '+')),
    'FBtr0079531_100_50': ('2L:7986605-7986755(+)', 0, ('2L', 7986705, '+')),
    'FBtr0079531_50_100': ('2L:7986655-7986785(+)', 0, ('2L', 7986705, '+')),
    'FBtr0079531_at_splice_0_100': ('2L:7986774-7986785(+)',0,('2L', 7986774, '+')),
    'FBtr0079531_at_splice_100_0': ('2L:7986605-7986705(+)',0,('2L', 7986774, '+')),
    'FBtr0079531_at_splice_100_50': ('2L:7986605-7986705^7986774-7986785(+)',0,('2L', 7986774, '+')),
    'FBtr0079531_at_splice_50_100': ('2L:7986655-7986705^7986774-7986785(+)',0,('2L', 7986774, '+')),
    'FBtr0079531_no_cds_0_100': ('na', nan, nan),
    'FBtr0079531_no_cds_100_0': ('na', nan, nan),
    'FBtr0079531_no_cds_100_50': ('na', nan, nan),
    'FBtr0079531_no_cds_50_100': ('na', nan, nan),
    'FBtr0079531_no_utr_0_100': ('na', 0, ('2L', 7986705, '+')),
    'FBtr0079531_no_utr_100_0': ('2L:7986605-7986705(+)',0,('2L', 7986705, '+')),
    'FBtr0079531_no_utr_100_50': ('2L:7986605-7986705(+)',0,('2L', 7986705, '+')),
    'FBtr0079531_no_utr_50_100': ('2L:7986655-7986705(+)',0,('2L', 7986705, '+')),
    'FBtr0079531_short_utr_0_100': ('2L:7986705-7986715(+)',0,('2L', 7986705, '+')),
    'FBtr0079531_short_utr_100_0': ('2L:7986605-7986705(+)',0,('2L', 7986705, '+')),
    'FBtr0079531_short_utr_100_50': ('2L:7986605-7986715(+)',0,('2L', 7986705, '+')),
    'FBtr0079531_short_utr_50_100': ('2L:7986655-7986715(+)',0,('2L', 7986705, '+')),
    'FBtr0081950_0_100': ('3R:4517499-4517599(-)', 0, ('3R', 4517598, '-')),
    'FBtr0081950_100_0': ('3R:4517599-4517699(-)', 0, ('3R', 4517598, '-')),
    'FBtr0081950_100_50': ('3R:4517549-4517699(-)', 0, ('3R', 4517598, '-')),
    'FBtr0081950_50_100': ('3R:4517499-4517649(-)', 0, ('3R', 4517598, '-')),
    'FBtr0081950_at_splice_0_100': ('3R:4517210-4517241(-)',0,('3R', 4517240, '-')),
    'FBtr0081950_at_splice_100_0': ('3R:4517599-4517699(-)',0,('3R', 4517240, '-')),
    'FBtr0081950_at_splice_100_50': ('3R:4517210-4517241^4517599-4517699(-)',0,('3R', 4517240, '-')),
    'FBtr0081950_at_splice_50_100': ('3R:4517210-4517241^4517599-4517649(-)',0,('3R', 4517240, '-')),
    'FBtr0081950_no_cds_0_100': ('na', nan, nan),
    'FBtr0081950_no_cds_100_0': ('na', nan, nan),
    'FBtr0081950_no_cds_100_50': ('na', nan, nan),
    'FBtr0081950_no_cds_50_100': ('na', nan, nan),
    'FBtr0081950_no_utr_0_100': ('na', 0, ('3R', 4517598, '-')),
    'FBtr0081950_no_utr_100_0': ('3R:4517599-4517699(-)',0,('3R', 4517598, '-')),
    'FBtr0081950_no_utr_100_50': ('3R:4517599-4517699(-)',0,('3R', 4517598, '-')),
    'FBtr0081950_no_utr_50_100': ('3R:4517599-4517649(-)',0,('3R', 4517598, '-')),
    'FBtr0081950_short_utr_0_100': ('3R:4517589-4517599(-)',0,('3R', 4517598, '-')),
    'FBtr0081950_short_utr_100_0': ('3R:4517599-4517699(-)',0,('3R', 4517598, '-')),
    'FBtr0081950_short_utr_100_50': ('3R:4517589-4517699(-)',0,('3R', 4517598, '-')),
    'FBtr0081950_short_utr_50_100': ('3R:4517589-4517649(-)',0,('3R', 4517598, '-'))
}

_DO_GENERATE_MAX_WINDOW = {
    "3_same_start_plus"      : ["FBtr0306336","FBtr0079531","FBtr0079531_short_utr"],
    "3_same_start_plus2"     : ["FBtr0306336","FBtr0079531","FBtr0079531_no_utr"],
    "3_diff_start_plus"      : ["FBtr0306336","FBtr0079531","FBtr0079531_short_utr","FBtr0300185"],
    "3_same_plus_ncrna_plus" : ["FBtr0306336","FBtr0079531","FBtr0079531_short_utr","FBtr0079531_no_cds"],
    "3_diff_plus_ncrna_plus" : ["FBtr0306336","FBtr0079531","FBtr0079531_short_utr","FBtr0300185","FBtr0079531_no_cds"],

    "3_same_start_minus"       : ["FBtr0081950","FBtr0081951","FBtr0081950_short_utr"],
    "3_same_start_minus2"      : ["FBtr0081950","FBtr0081951","FBtr0081950_no_utr"],
    "3_diff_start_minus"       : ["FBtr0081950","FBtr0081951_alt_start","FBtr0081951"],
    "3_same_minus_ncrna_minus" : ["FBtr0081950","FBtr0081951","FBtr0081950_short_utr","FBtr0081950_no_cds"],
    "3_diff_minus_ncrna_minus" : ["FBtr0081950","FBtr0081951_alt_start","FBtr0081951","FBtr0081950_no_cds"],


    "stop_max_spanning_window_threeprime_plus"  : ["FBtr0079763","FBtr0303880"],

    "stop_max_spanning_window_threeprime_minus" : ["FBtr0079813","FBtr0303900"],

}

_DO_GENERATE_MAX_WINDOW_RESULTS = {
    '3_diff_minus_ncrna_minus_0_100': ['na', nan, nan],
    '3_diff_minus_ncrna_minus_100_0': ['na', nan, nan],
    '3_diff_minus_ncrna_minus_100_50': ['na', nan, nan],
    '3_diff_minus_ncrna_minus_50_100': ['na', nan, nan],
    '3_diff_plus_ncrna_plus_0_100': ['na', nan, nan],
    '3_diff_plus_ncrna_plus_100_0': ['na', nan, nan],
    '3_diff_plus_ncrna_plus_100_50': ['na', nan, nan],
    '3_diff_plus_ncrna_plus_50_100': ['na', nan, nan],
    '3_diff_start_minus_0_100': ['na', nan, nan],
    '3_diff_start_minus_100_0': ['na', nan, nan],
    '3_diff_start_minus_100_50': ['na', nan, nan],
    '3_diff_start_minus_50_100': ['na', nan, nan],
    '3_diff_start_plus_0_100': ['na', nan, nan],
    '3_diff_start_plus_100_0': ['na', nan, nan],
    '3_diff_start_plus_100_50': ['na', nan, nan],
    '3_diff_start_plus_50_100': ['na', nan, nan],
    '3_same_minus_ncrna_minus_0_100': ['na', nan, nan],
    '3_same_minus_ncrna_minus_100_0': ['na', nan, nan],
    '3_same_minus_ncrna_minus_100_50': ['na', nan, nan],
    '3_same_minus_ncrna_minus_50_100': ['na', nan, nan],
    '3_same_plus_ncrna_plus_0_100': ['na', nan, nan],
    '3_same_plus_ncrna_plus_100_0': ['na', nan, nan],
    '3_same_plus_ncrna_plus_100_50': ['na', nan, nan],
    '3_same_plus_ncrna_plus_50_100': ['na', nan, nan],
    '3_same_start_minus2_0_100': ['3R:4519776-4519876(-)', 0, 0],
    '3_same_start_minus2_100_0': ['na', nan, nan],
    '3_same_start_minus2_100_50': ['3R:4519826-4519876(-)', 100, 100],
    '3_same_start_minus2_50_100': ['3R:4519776-4519876(-)', 50, 50],
    '3_same_start_minus_0_100': ['3R:4519776-4519876(-)', 0, 0],
    '3_same_start_minus_100_0': ['3R:4519876-4519894(-)', 82, 100],
    '3_same_start_minus_100_50': ['3R:4519826-4519894(-)', 82, 100],
    '3_same_start_minus_50_100': ['3R:4519776-4519894(-)', 32, 50],
    '3_same_start_plus2_0_100': ['2L:7985674-7985768^7985833-7985839(+)', 0, 0],
    '3_same_start_plus2_100_0': ['na', nan, nan],
    '3_same_start_plus2_100_50': ['2L:7985674-7985724(+)', 100, 100],
    '3_same_start_plus2_50_100': ['2L:7985674-7985768^7985833-7985839(+)',50,50],
    '3_same_start_plus_0_100': ['2L:7985674-7985768^7985833-7985839(+)', 0, 0],
    '3_same_start_plus_100_0': ['2L:7985664-7985674(+)', 90, 100],
    '3_same_start_plus_100_50': ['2L:7985664-7985724(+)', 90, 100],
    '3_same_start_plus_50_100': ['2L:7985664-7985768^7985833-7985839(+)', 40, 50],
    
    'stop_max_spanning_window_threeprime_minus_0_100': ['2L:9397646-9397746(-)',  0,   0],
    'stop_max_spanning_window_threeprime_minus_100_0': ['2L:9397746-9397764(-)',  0,   100],
    'stop_max_spanning_window_threeprime_minus_100_50': ['2L:9397696-9397764(-)', 82,  100],
    'stop_max_spanning_window_threeprime_minus_50_100': ['2L:9397646-9397764(-)', 32,  50],
    
    'stop_max_spanning_window_threeprime_plus_0_100': ['2L:8997663-8997668(+)',   0,   0],
    'stop_max_spanning_window_threeprime_plus_100_0': ['2L:8997640-8997663(+)',   77,  100],
    'stop_max_spanning_window_threeprime_plus_100_50': ['2L:8997640-8997668(+)',  77,  100],
    'stop_max_spanning_window_threeprime_plus_50_100': ['2L:8997640-8997668(+)',  27,  50],  
}

_DO_GENERATE_MAX_WINDOW_RESULTS_MASKED = {
    '3_diff_minus_ncrna_minus_0_100': ['na', nan, nan, 'na'],
    '3_diff_minus_ncrna_minus_100_0': ['na', nan, nan, 'na'],
    '3_diff_minus_ncrna_minus_100_50': ['na', nan, nan, 'na'],
    '3_diff_minus_ncrna_minus_50_100': ['na', nan, nan, 'na'],
    '3_diff_plus_ncrna_plus_0_100': ['na', nan, nan, 'na'],
    '3_diff_plus_ncrna_plus_100_0': ['na', nan, nan, 'na'],
    '3_diff_plus_ncrna_plus_100_50': ['na', nan, nan, 'na'],
    '3_diff_plus_ncrna_plus_50_100': ['na', nan, nan, 'na'],
    '3_diff_start_minus_0_100': ['na', nan, nan, 'na'],
    '3_diff_start_minus_100_0': ['na', nan, nan, 'na'],
    '3_diff_start_minus_100_50': ['na', nan, nan, 'na'],
    '3_diff_start_minus_50_100': ['na', nan, nan, 'na'],
    '3_diff_start_plus_0_100': ['na', nan, nan, 'na'],
    '3_diff_start_plus_100_0': ['na', nan, nan, 'na'],
    '3_diff_start_plus_100_50': ['na', nan, nan, 'na'],
    '3_diff_start_plus_50_100': ['na', nan, nan, 'na'],
    '3_same_minus_ncrna_minus_0_100': ['na', nan, nan, 'na'],
    '3_same_minus_ncrna_minus_100_0': ['na', nan, nan, 'na'],
    '3_same_minus_ncrna_minus_100_50': ['na', nan, nan, 'na'],
    '3_same_minus_ncrna_minus_50_100': ['na', nan, nan, 'na'],
    '3_same_plus_ncrna_plus_0_100': ['na', nan, nan, 'na'],
    '3_same_plus_ncrna_plus_100_0': ['na', nan, nan, 'na'],
    '3_same_plus_ncrna_plus_100_50': ['na', nan, nan, 'na'],
    '3_same_plus_ncrna_plus_50_100': ['na', nan, nan, 'na'],
    '3_same_start_minus2_0_100': ['3R:4519776-4519876(-)', 0, 0, 'na'],
    '3_same_start_minus2_100_0': ['na', nan, nan, 'na'],
    '3_same_start_minus2_100_50': ['3R:4519826-4519876(-)', 100, 100, 'na'],
    '3_same_start_minus2_50_100': ['3R:4519776-4519876(-)', 50, 50, 'na'],
    '3_same_start_minus_0_100': ['3R:4519776-4519876(-)', 0, 0, 'na'],
    '3_same_start_minus_100_0': ['3R:4519876-4519894(-)',82,100,'3R:4519879-4519891(-)'],
    '3_same_start_minus_100_50': ['3R:4519826-4519894(-)',82,100,'3R:4519879-4519891(-)'],
    '3_same_start_minus_50_100': ['3R:4519776-4519894(-)',32,50,'3R:4519879-4519891(-)'],
    '3_same_start_plus2_0_100': ['2L:7985674-7985768^7985833-7985839(+)',0,0,'2L:7985694-7985744(+)'],
    '3_same_start_plus2_100_0': ['na', nan, nan, 'na'],
    '3_same_start_plus2_100_50': ['2L:7985674-7985724(+)',100,100,'2L:7985694-7985724(+)'],
    '3_same_start_plus2_50_100': ['2L:7985674-7985768^7985833-7985839(+)',50,50,'2L:7985694-7985744(+)'],
    '3_same_start_plus_0_100': ['2L:7985674-7985768^7985833-7985839(+)',0,0,'2L:7985694-7985744(+)'],
    '3_same_start_plus_100_0': ['2L:7985664-7985674(+)', 90, 100, 'na'],
    '3_same_start_plus_100_50': ['2L:7985664-7985724(+)',90,100,'2L:7985694-7985724(+)'],
    '3_same_start_plus_50_100': ['2L:7985664-7985768^7985833-7985839(+)',40,50,'2L:7985694-7985744(+)'],
    
    
    'stop_max_spanning_window_threeprime_minus_0_100': ['2L:9397646-9397746(-)',  0,   0],
    'stop_max_spanning_window_threeprime_minus_100_0': ['2L:9397746-9397764(-)',  0,   100],
    'stop_max_spanning_window_threeprime_minus_100_50': ['2L:9397696-9397764(-)', 82,  100],
    'stop_max_spanning_window_threeprime_minus_50_100': ['2L:9397646-9397764(-)', 32,  50],
    
    'stop_max_spanning_window_threeprime_plus_0_100': ['2L:8997663-8997668(+)',   0,   0],
    'stop_max_spanning_window_threeprime_plus_100_0': ['2L:8997640-8997663(+)',   77,  100],
    'stop_max_spanning_window_threeprime_plus_100_50': ['2L:8997640-8997668(+)',  77,  100],
    'stop_max_spanning_window_threeprime_plus_50_100': ['2L:8997640-8997668(+)',  27,  50],      
}

_DO_GENERATE_MULTI_GENE = {
    "3_same_start"         : _DO_GENERATE_MAX_WINDOW["3_same_start_plus"] + _DO_GENERATE_MAX_WINDOW["3_same_start_minus"],
    "3_same_start_1_ncrna" : _DO_GENERATE_MAX_WINDOW["3_same_plus_ncrna_plus"] + _DO_GENERATE_MAX_WINDOW["3_same_start_minus"],
    "same_and_diff"        : _DO_GENERATE_MAX_WINDOW["3_diff_start_plus"] + _DO_GENERATE_MAX_WINDOW["3_same_start_minus"]
}

_DO_GENERATE_MULTI_GENE_RESULTS = {
    "3_same_start_50_100"         : [_DO_GENERATE_MAX_WINDOW_RESULTS["3_same_start_minus_50_100"],
                                     _DO_GENERATE_MAX_WINDOW_RESULTS["3_same_start_plus_50_100"],
                                     ],
    "3_same_start_1_ncrna_50_100" : [_DO_GENERATE_MAX_WINDOW_RESULTS["3_same_plus_ncrna_plus_50_100"],
                                     _DO_GENERATE_MAX_WINDOW_RESULTS["3_same_start_minus_50_100"]],
    "same_and_diff_50_100"        : [_DO_GENERATE_MAX_WINDOW_RESULTS["3_diff_start_plus_50_100"],
                                     _DO_GENERATE_MAX_WINDOW_RESULTS["3_same_start_minus_50_100"]],
}
