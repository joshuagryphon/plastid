#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.readers.bed`

    :py:func:`BED_to_Transcripts`
        Read BED files to Transcript objects
    
    :py:func:`BED_to_SegmentChain`
        Reads BED files to SegmentChain objects
    
    :py:class:`BED_Reader`
        Reads BED files to SegmentChain objects

See http://genome.ucsc.edu/FAQ/FAQformat.html
"""
import unittest
import functools
import warnings
import sys
import pandas as pd

from csv import QUOTE_NONE
from nose.plugins.attrib import attr
from plastid.util.services.mini2to3 import cStringIO
from plastid.genomics.roitools import SegmentChain, GenomicSegment, Transcript
from plastid.readers.bed import BED_to_Transcripts, BED_to_SegmentChain, BED_Reader
from nose.tools import assert_equal, assert_true, assert_dict_equal, assert_greater_equal

warnings.simplefilter("ignore",DeprecationWarning)

#===============================================================================
# INDEX: test suites
#===============================================================================
@attr(test="unit")
class TestBED():
    """Test case for BED input/output"""
    
    @classmethod
    def setUpClass(cls):
        cls.header = _BED_HEADER
        cls.data = {}
        cls.extracol_data = {}
        bed_df = pd.read_table(cStringIO.StringIO(_BED12_DATA),header=None,sep="\t",index_col=None)
        extra_df = pd.read_table(cStringIO.StringIO(_EXTRA_COLS),header=0,sep="\t",index_col=None)
        cls.big_df = pd.concat([bed_df,extra_df],axis=1)

        for n in (3,4,5,6,8,9,12):
            cls.data[n] = cls.get_bed_subset(cls.header,n,0)
            cls.extracol_data[n] = cls.get_bed_subset(cls.header,n,4)

    @classmethod
    def get_bed_subset(cls,header,bed_cols,extra_cols=0):
        buf = cStringIO.StringIO()
        columns = cls.big_df.columns[list(range(bed_cols)) + list(range(12,12+extra_cols))]
        cls.big_df.to_csv(buf,columns=columns,sep="\t",index=False,header=False,quoting=QUOTE_NONE) #,float_format="%.8f")
        return buf.getvalue()
   
    @staticmethod
    def check_equal(found,expected,msg=None):
        if msg is not None:
            assert_equal(found,expected,msg)
        else:
            assert_equal(found,expected)

    def test_bed_import_3to12plus4_columns_with_formatters(self):
        names = [("numcol",int),
                 ("floatcol",float),
                 ("strcol",str),
                 ("attrcol",str),
                ]

        tx_reader = functools.partial(BED_Reader,return_type=Transcript,extra_columns=names)
        seg_reader = functools.partial(BED_Reader,return_type=SegmentChain,extra_columns=names)
        tests = [(seg_reader,_TEST_SEGMENTCHAINS,"reader_segmentchain"),
                 (tx_reader,_TEST_TRANSCRIPTS,"reader_transcript"),
                ]
        for reader_fn, known_set, name in tests:
            for n,data_str in sorted(self.extracol_data.items()):
                c = 0
                for (test_ivc,known_ivc) in zip(reader_fn(cStringIO.StringIO(data_str)),
                                                           known_set):
                    for x in range(4):
                        colname = names[x][0]
                        assert_true(colname in test_ivc.attr,"Column name '%s' not found in attr dict (%s BED columns)" % (x,n))
                        assert_equal(test_ivc.attr[colname],self.big_df.iloc[c,12+x])

                    # columns: chrom, start, end
                    if n >= 3:
                        # no strand info, so we need to test iv.start, iv.end, iv.chrom
                        err_msg = "%s failed endpoint equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.spanning_segment.start,test_ivc.spanning_segment.start,err_msg
                        yield self.check_equal, known_ivc.spanning_segment.end,test_ivc.spanning_segment.end,err_msg
                        yield self.check_equal, known_ivc.spanning_segment.chrom,test_ivc.spanning_segment.chrom,err_msg
                    # column: name
                    if n >= 4:
                        err_msg = "%s failed name equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr["ID"],test_ivc.attr["ID"],err_msg
                    # column: score
                    if n >= 5:
                        err_msg = "%s failed score equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr.get("score",0),test_ivc.attr["score"],err_msg
                    # column : strand
                    if n >= 6:
                        err_msg = "%s failed strand equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.spanning_segment.strand,test_ivc.spanning_segment.strand
                    # column: color
                    if n >= 9:
                        err_msg = "%s failed color equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr.get("color","#000000"),test_ivc.attr["color"],err_msg
                    # columns: exon/block info
                    if n == 12:
                        err_msg = "%s failed block equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        for iv1, iv2 in zip(known_ivc,test_ivc):
                            assert_equal(iv1,iv2,err_msg)
                        err_msg = "%s failed position set on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.get_position_set(), test_ivc.get_position_set(), err_msg
                    
                    c += 1
                
                yield self.check_equal, c,len(known_set),"Not all intervals loaded! Expected %s, found %s." % (len(known_set),c)

    def test_bed_import_3to12plus4_columns_with_names(self):
        names = [X for X in self.big_df.columns[-4:]]
        tx_reader = functools.partial(BED_Reader,return_type=Transcript,extra_columns=names)
        seg_reader = functools.partial(BED_Reader,return_type=SegmentChain,extra_columns=names)
        tests = [(seg_reader,_TEST_SEGMENTCHAINS,"reader_segmentchain"),
                 (tx_reader,_TEST_TRANSCRIPTS,"reader_transcript"),
                ]
        for reader_fn, known_set, name in tests:
            for n,data_str in sorted(self.extracol_data.items()):
                c = 0
                for (test_ivc,known_ivc) in zip(reader_fn(cStringIO.StringIO(data_str)),
                                                           known_set):
                    for x in range(4):
                        colname = names[x]
                        assert_true(colname in test_ivc.attr)
                        assert_equal(str(test_ivc.attr[colname]),str(self.big_df.iloc[c,12+x]))

                    # columns: chrom, start, end
                    if n >= 3:
                        # no strand info, so we need to test iv.start, iv.end, iv.chrom
                        err_msg = "%s failed endpoint equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.spanning_segment.start,test_ivc.spanning_segment.start,err_msg
                        yield self.check_equal, known_ivc.spanning_segment.end,test_ivc.spanning_segment.end,err_msg
                        yield self.check_equal, known_ivc.spanning_segment.chrom,test_ivc.spanning_segment.chrom,err_msg
                    # column: name
                    if n >= 4:
                        err_msg = "%s failed name equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr["ID"],test_ivc.attr["ID"],err_msg
                    # column: score
                    if n >= 5:
                        err_msg = "%s failed score equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr.get("score",0),test_ivc.attr["score"],err_msg
                    # column : strand
                    if n >= 6:
                        err_msg = "%s failed strand equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.spanning_segment.strand,test_ivc.spanning_segment.strand
                    # column: color
                    if n >= 9:
                        err_msg = "%s failed color equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr.get("color","#000000"),test_ivc.attr["color"],err_msg
                    # columns: exon/block info
                    if n == 12:
                        err_msg = "%s failed block equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        for iv1, iv2 in zip(known_ivc,test_ivc):
                            assert_equal(iv1,iv2,err_msg)
                        err_msg = "%s failed position set on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.get_position_set(), test_ivc.get_position_set(), err_msg
                    
                    c += 1
                
                yield self.check_equal, c,len(known_set),"Not all intervals loaded! Expected %s, found %s." % (len(known_set),c)

    def test_bed_import_3to12plus4_columns_with_int(self):
        tx_reader = functools.partial(BED_Reader,return_type=Transcript,extra_columns=4)
        seg_reader = functools.partial(BED_Reader,return_type=SegmentChain,extra_columns=4)
        tests = [(seg_reader,_TEST_SEGMENTCHAINS,"reader_segmentchain"),
                 (tx_reader,_TEST_TRANSCRIPTS,"reader_transcript"),
                ]
        for reader_fn, known_set, name in tests:
            for n,data_str in sorted(self.extracol_data.items()):
                c = 0
                for (test_ivc,known_ivc) in zip(reader_fn(cStringIO.StringIO(data_str)),
                                                           known_set):
                    for x in range(4):
                        colname = "custom%s" % x
                        assert_true(colname in test_ivc.attr)
                        assert_equal(str(test_ivc.attr[colname]),str(self.big_df.iloc[c,12+x]))

                    # columns: chrom, start, end
                    if n >= 3:
                        # no strand info, so we need to test iv.start, iv.end, iv.chrom
                        err_msg = "%s failed endpoint equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.spanning_segment.start,test_ivc.spanning_segment.start,err_msg
                        yield self.check_equal, known_ivc.spanning_segment.end,test_ivc.spanning_segment.end,err_msg
                        yield self.check_equal, known_ivc.spanning_segment.chrom,test_ivc.spanning_segment.chrom,err_msg
                    # column: name
                    if n >= 4:
                        err_msg = "%s failed name equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr["ID"],test_ivc.attr["ID"],err_msg
                    # column: score
                    if n >= 5:
                        err_msg = "%s failed score equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr.get("score",0),test_ivc.attr["score"],err_msg
                    # column : strand
                    if n >= 6:
                        err_msg = "%s failed strand equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.spanning_segment.strand,test_ivc.spanning_segment.strand
                    # column: color
                    if n >= 9:
                        err_msg = "%s failed color equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr.get("color","#000000"),test_ivc.attr["color"],err_msg
                    # columns: exon/block info
                    if n == 12:
                        err_msg = "%s failed block equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        for iv1, iv2 in zip(known_ivc,test_ivc):
                            assert_equal(iv1,iv2,err_msg)
                        err_msg = "%s failed position set on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.get_position_set(), test_ivc.get_position_set(), err_msg
                    
                    c += 1
                
                yield self.check_equal, c,len(known_set),"Not all intervals loaded! Expected %s, found %s." % (len(known_set),c)

    def test_bed_export_3to12plus4_columns_with_names(self):
        names = [X for X in self.big_df.columns[-4:]]
        tests = [(Transcript.from_bed,_TEST_TRANSCRIPTS,"tx_frombed_plus4_int"),
                 (SegmentChain.from_bed,_TEST_SEGMENTCHAINS,"segchain_frombed_plus4_int"),
                ]
        for import_fn, known_set, name in tests:
            extracol12plus = [X.split("\t") for X in self.extracol_data[12].strip("\n").split("\n")]
            for n,data_str in sorted(self.extracol_data.items()):
                for c, line in enumerate(data_str.strip("\n").split("\n")):
                    out_line = import_fn(line,extra_columns=names).as_bed(as_int=False).strip("\n")
                    out_items = out_line.split("\t")[:n] + out_line.split("\t")[-4:]
                    expected_items = extracol12plus[c][:n] + extracol12plus[c][-4:]
                    msg = "%s BED %s+%s export unequal for lines:\nin:  %s\nout: %s\nexp: %s" % (name, n,4,line,
                            "\t".join(out_items),"\t".join(expected_items))
                    yield self.check_equal, out_items, expected_items, msg
                
                yield self.check_equal, c+1, len(known_set),"Not all intervals loaded! Expected %s, found %s." % (len(known_set),c)

    def test_bed_export_3to12plus4_columns_with_int(self):
        tests = [(Transcript.from_bed,_TEST_TRANSCRIPTS,"tx_frombed_plus4_int"),
                 (SegmentChain.from_bed,_TEST_SEGMENTCHAINS,"segchain_frombed_plus4_int"),
                ]
        for import_fn, known_set, name in tests:
            extracol12plus = [X.split("\t") for X in self.extracol_data[12].strip("\n").split("\n")]
            for n,data_str in sorted(self.extracol_data.items()):
                for c, line in enumerate(data_str.strip("\n").split("\n")):
                    out_line = import_fn(line,extra_columns=4).as_bed(as_int=False).strip("\n")
                    out_items = out_line.split("\t")[:n] + out_line.split("\t")[-4:]
                    expected_items = extracol12plus[c][:n] + extracol12plus[c][-4:]
                    msg = "%s BED %s+%s export unequal for lines:\nin:  %s\nout: %s\nexp: %s" % (name, n,4,line,
                            "\t".join(out_items),"\t".join(expected_items))
                    yield self.check_equal, out_items, expected_items, msg
                
                yield self.check_equal, c+1, len(known_set),"Not all intervals loaded! Expected %s, found %s." % (len(known_set),c)

    def test_bed_import_3to12_columns(self):
        tx_reader = functools.partial(BED_Reader,return_type=Transcript)
        tests = [(BED_to_SegmentChain,_TEST_SEGMENTCHAINS,"tosegmentchain,segmentchain"),
                 (BED_to_Transcripts,_TEST_TRANSCRIPTS,"totranscripts_transcript"),
                 (BED_Reader,_TEST_SEGMENTCHAINS,"reader_segmentchain"),
                 (tx_reader,_TEST_TRANSCRIPTS,"reader_transcript"),
                ]
        for reader_fn, known_set, name in tests:
            for n,data_str in sorted(self.data.items()):
                c = 0
                for (test_ivc,known_ivc) in zip(reader_fn(cStringIO.StringIO(data_str)),
                                                           known_set):
                    # columns: chrom, start, end
                    if n >= 3:
                        # no strand info, so we need to test iv.start, iv.end, iv.chrom
                        err_msg = "%s failed endpoint equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.spanning_segment.start,test_ivc.spanning_segment.start,err_msg
                        yield self.check_equal, known_ivc.spanning_segment.end,test_ivc.spanning_segment.end,err_msg
                        yield self.check_equal, known_ivc.spanning_segment.chrom,test_ivc.spanning_segment.chrom,err_msg
                    # column: name
                    if n >= 4:
                        err_msg = "%s failed name equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr["ID"],test_ivc.attr["ID"],err_msg
                    # column: score
                    if n >= 5:
                        err_msg = "%s failed score equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr.get("score",0),test_ivc.attr["score"],err_msg
                    # column : strand
                    if n >= 6:
                        err_msg = "%s failed strand equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.spanning_segment.strand,test_ivc.spanning_segment.strand
                    # column: color
                    if n >= 9:
                        err_msg = "%s failed color equality on %s-column BED input: %s,%s" % (name, n,known_ivc.attr,test_ivc.attr)
                        yield self.check_equal, known_ivc.attr.get("color","#000000"),test_ivc.attr["color"],err_msg
                    # columns: exon/block info
                    if n == 12:
                        err_msg = "%s failed block equality on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        for iv1, iv2 in zip(known_ivc,test_ivc):
                            assert_equal(iv1,iv2,err_msg)
                        err_msg = "%s failed position set on %s-column BED input: %s,%s" % (name, n,known_ivc,test_ivc)
                        yield self.check_equal, known_ivc.get_position_set(), test_ivc.get_position_set(), err_msg
                    
                    c += 1
                
                yield self.check_equal, c,len(known_set),"Not all intervals loaded! Expected %s, found %s." % (len(known_set),c)

    def test_ivcollection_thick_start_end_8to12_columns(self):
        """Checks equality of thickstart and thickend attributes for SegmentChain objects"""
        for n,data_str in sorted(self.data.items()):
            for c, (test_ivc,known_ivc) in enumerate(zip(BED_Reader(cStringIO.StringIO(data_str),return_type=SegmentChain),
                                                       _TEST_SEGMENTCHAINS)):
                if n >= 8:
                    err_msg = "Failed thickstart/end equality on %s-column BED input: %s,%s" % (n,known_ivc.attr,test_ivc.attr)
                    if known_ivc.attr.get("thickstart",None) is not None:
                        yield self.check_equal, known_ivc.attr["thickstart"],test_ivc.attr["thickstart"],err_msg
                    if known_ivc.attr.get("thickend",None) is not None:
                        yield self.check_equal, known_ivc.attr.get("thickend"),test_ivc.attr["thickend"],err_msg
        
            yield self.check_equal, c,20-1,"Not all intervals loaded! Expected %s, found %s." % (20-1,c)

    def test_transcript_cds_start_end_8to12_columns(self):
        """Checks equality of endpoints of coding regions for Transcript objects"""
        for n,data_str in sorted(self.data.items()):
            for c, (test_ivc,known_ivc) in enumerate(zip(BED_Reader(cStringIO.StringIO(data_str),return_type=Transcript),
                                                       _TEST_TRANSCRIPTS)):
                if n >= 8:
                    err_msg = "Failed thickstart/end equality on %s-column BED input: %s,%s" % (n,known_ivc.attr,test_ivc.attr)
                    if known_ivc.attr.get("cds_genome_start",None) is not None:
                        yield self.check_equal, known_ivc.attr["cds_start"],test_ivc.attr["cds_start"],err_msg
                        yield self.check_equal, known_ivc.attr["cds_genome_start"],test_ivc.attr["cds_genome_start"],err_msg
                        yield self.check_equal, known_ivc.cds_genome_start,test_ivc.cds_genome_start,err_msg
                        yield self.check_equal, known_ivc.cds_start,test_ivc.cds_start,err_msg
                    if known_ivc.attr.get("cds_genome_end",None) is not None:
                        yield self.check_equal, known_ivc.attr["cds_end"],test_ivc.attr["cds_end"],err_msg
                        yield self.check_equal, known_ivc.attr["cds_genome_end"],test_ivc.attr["cds_genome_end"],err_msg
                        yield self.check_equal, known_ivc.cds_genome_end,test_ivc.cds_genome_end,err_msg
                        yield self.check_equal, known_ivc.cds_end,test_ivc.cds_end,err_msg
            
            yield self.check_equal, c,20-1,"Not all intervals loaded! Expected %s, found %s." % (20-1,c)
        
    def test_track_subtype_parsing(self):
        reader = BED_Reader(cStringIO.StringIO(_NARROW_PEAK_TEXT))
        for c, (found, expected) in enumerate(zip(reader,_NARROW_PEAK_CHAINS)):
            found.attr.pop("color")
            found.attr.pop("score")
            assert_dict_equal(found.attr,expected.attr)
            assert_equal(found,expected)

        assert_equal(c,len(_NARROW_PEAK_CHAINS)-1)

    def test_track_subtype_raises_warning_if_wrong_extra_columns(self):
        reader = BED_Reader(cStringIO.StringIO(_NARROW_PEAK_TEXT),extra_columns=14)
        with warnings.catch_warnings(record=True) as warns:
            warnings.simplefilter("always")
            ltmp = list(reader)
            assert_greater_equal(len(warns),0)

#===============================================================================
# INDEX: test data
#===============================================================================

# test dataset, constructed manually to include various edge cases
_TEST_SEGMENTCHAINS = [
    # single-interval
    SegmentChain(GenomicSegment("chrA",100,1100,"+"),ID="IVC1p"),
    SegmentChain(GenomicSegment("chrA",100,1100,"-"),ID="IVC1m"),
    # multi-interval
    SegmentChain(GenomicSegment("chrA",100,1100,"+"),GenomicSegment("chrA",2100,2600,"+"),ID="IVC2p"),
    SegmentChain(GenomicSegment("chrA",100,1100,"-"),GenomicSegment("chrA",2100,2600,"-"),ID="IVC2m"),
    # multi-interval, with score
    SegmentChain(GenomicSegment("chrA",100,1100,"+"),GenomicSegment("chrA",2100,2600,"+"),ID="IVC3p",score=500),
    SegmentChain(GenomicSegment("chrA",100,1100,"-"),GenomicSegment("chrA",2100,2600,"-"),ID="IVC3m",score=500),
    SegmentChain(GenomicSegment("chrA",100,1100,"+"),GenomicSegment("chrA",2100,2600,"+"),GenomicSegment("chrA",2605,2700,"+"),
                 ID="IVC4p",score=500),
    SegmentChain(GenomicSegment("chrA",100,1100,"-"),GenomicSegment("chrA",2100,2600,"-"),
                 ID="IVC4m",score=500),
    # multi-interval, with score and color
    SegmentChain(GenomicSegment("chrA",100,1100,"+"),GenomicSegment("chrA",2100,2600,"+"),GenomicSegment("chrA",2605,2700,"+"),
                 ID="IVC5p",score=500,color="#007ADF"),
    SegmentChain(GenomicSegment("chrA",100,1100,"-"),GenomicSegment("chrA",2100,2600,"-"),GenomicSegment("chrA",2605,2700,"-"),
                 ID="IVC5m",score=500,color="#007ADF"),
    SegmentChain(GenomicSegment("chrA",100,1100,"+"),GenomicSegment("chrA",2100,2600,"+"),GenomicSegment("chrA",2605,2700,"+"),
                 ID="IVC6p",score=500,color="#007ADF"),
    SegmentChain(GenomicSegment("chrA",100,1100,"-"),GenomicSegment("chrA",2100,2600,"-"),GenomicSegment("chrA",2605,2700,"-"),
                 ID="IVC6m",score=500,color="#007ADF"),
    # multi-interval, with score, color, thickstart, and thickend, internally
    SegmentChain(GenomicSegment("chrA",100,1100,"+"),GenomicSegment("chrA",2100,2600,"+"),GenomicSegment("chrA",2605,2700,"+"),
                 ID="IVC7p",score=500,color="#007ADF",thickstart=2200,thickend=2400),
    SegmentChain(GenomicSegment("chrA",100,1100,"-"),GenomicSegment("chrA",2100,2600,"-"),GenomicSegment("chrA",2605,2700,"-"),
                 ID="IVC7m",score=500,color="#007ADF",thickstart=2200,thickend=2400),
    # multi-interval, thickend and thickstart covering whole SegmentChain             
    SegmentChain(GenomicSegment("chrA",100,1100,"+"),GenomicSegment("chrA",2100,2600,"+"),GenomicSegment("chrA",2605,2700,"+"),
                 ID="IVC8p",score=500,color="#007ADF",thickstart=100,thickend=2700),
    SegmentChain(GenomicSegment("chrA",100,1100,"-"),GenomicSegment("chrA",2100,2600,"-"),GenomicSegment("chrA",2605,2700,"-"),
                 ID="IVC8m",score=500,color="#007ADF",thickstart=100,thickend=2700),
    # multi-interval, thickend and thickstart at exon-exon junctions             
    SegmentChain(GenomicSegment("chrA",100,1100,"+"),GenomicSegment("chrA",2100,2600,"+"),GenomicSegment("chrA",2605,2700,"+"),
                 ID="IVC9p",score=500,color="#007ADF",thickstart=2100,thickend=2600),
    SegmentChain(GenomicSegment("chrA",100,1100,"-"),GenomicSegment("chrA",2100,2600,"-"),GenomicSegment("chrA",2605,2700,"-"),
                 ID="IVC9m",score=500,color="#007ADF",thickstart=2100,thickend=2600),
    # multi-interval, thickend and thickstart at exon-exon junctions             
    SegmentChain(GenomicSegment("chrA",100,1100,"+"),GenomicSegment("chrA",2100,2600,"+"),GenomicSegment("chrA",2605,2700,"+"),
                 ID="IVC10p",score=500,color="#007ADF",thickstart=1099,thickend=2101),
    SegmentChain(GenomicSegment("chrA",100,1100,"-"),GenomicSegment("chrA",2100,2600,"-"),GenomicSegment("chrA",2605,2700,"-"),
                 ID="IVC10m",score=500,color="#007ADF",thickstart=1099,thickend=2101),
]

# same data, as transcripts
_TEST_TRANSCRIPTS = [Transcript(*X._segments,**X.attr) for X in _TEST_SEGMENTCHAINS]


_BED_HEADER = """browser position chrA:100-1100
track name=test_data description='my test data'
"""

# same data, as BED12 block
_BED12_DATA = """chrA    100    1100    IVC1p    0.0    +    100    100    0,0,0    1    1000,    0,
chrA    100    1100    IVC1m    0.0    -    100    100    0,0,0    1    1000,    0,
chrA    100    2600    IVC2p    0.0    +    100    100    0,0,0    2    1000,500,    0,2000,
chrA    100    2600    IVC2m    0.0    -    100    100    0,0,0    2    1000,500,    0,2000,
chrA    100    2600    IVC3p    500.0    +    100    100    0,0,0    2    1000,500,    0,2000,
chrA    100    2600    IVC3m    500.0    -    100    100    0,0,0    2    1000,500,    0,2000,
chrA    100    2700    IVC4p    500.0    +    100    100    0,0,0    3    1000,500,95,    0,2000,2505,
chrA    100    2600    IVC4m    500.0    -    100    100    0,0,0    2    1000,500,    0,2000,
chrA    100    2700    IVC5p    500.0    +    100    100    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC5m    500.0    -    100    100    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC6p    500.0    +    100    100    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC6m    500.0    -    100    100    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC7p    500.0    +    2200    2400    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC7m    500.0    -    2200    2400    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC8p    500.0    +    100    2700    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC8m    500.0    -    100    2700    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC9p    500.0    +    2100    2600    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC9m    500.0    -    2100    2600    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC10p    500.0    +    1099    2101    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC10m    500.0    -    1099    2101    0,122,223    3    1000,500,95,    0,2000,2505,""".replace("    ","\t")

_EXTRA_COLS="""numcol    floatcol    strcol    attrcol
0    3.14    a    gene_id "gene_0"; transcript_id "transcript_0";
1    2.72523    abc    gene_id "gene_1"; transcript_id "transcript_1";
2    30.12350    DEF    gene_id "gene_2"; transcript_id "transcript_2";
3    15123.20    ghi    gene_id "gene_3"; transcript_id "transcript_3";
4    2.0    alongword    gene_id "gene_4"; transcript_id "transcript_4";
5    -3.1234    a sentence with spaces    gene_id "gene_5"; transcript_id "transcript_5";
6    -20.5    some notes with "quotes"    gene_id "gene_6"; transcript_id "transcript_6";
7    -1e10    1    gene_id "gene_7"; transcript_id "transcript_7";
8    2e5    2    gene_id "gene_8"; transcript_id "transcript_8";
9    2.3e6    3.0    gene_id "gene_9"; transcript_id "transcript_9";
10    0.03    string1    gene_id "gene_10"; transcript_id "transcript_10";
11    1.0    string2    gene_id "gene_11"; transcript_id "transcript_11";
12    2.0    string3    gene_id "gene_12"; transcript_id "transcript_12";
13    3.0    string4 string5 string6    gene_id "gene_13"; transcript_id "transcript_13";
14    4.0    test    gene_id "gene_14"; transcript_id "transcript_14";
15    5.0    testetst    gene_id "gene_15"; transcript_id "transcript_15";
16    6.0    testsatsdfasf    gene_id "gene_16"; transcript_id "transcript_16";
17    7.0    asdgahghfzgdasdfasdf    gene_id "gene_17"; transcript_id "transcript_17";
18    8.0    asdfasdfadsfgaasdg    gene_id "gene_18"; transcript_id "transcript_18";
19    9.0    asdfasdfdasfdas    gene_id "gene_19"; transcript_id "transcript_19";
""".replace("    ","\t")


_NARROW_PEAK_TEXT = """track type=narrowPeak
chrI    100    15000    feature1    0    +    341.2    -123.2    -513.3    50
chrII    320    15000    feature2    0    -    2.1    -5123.2    0    650""".replace("    ","\t")

_NARROW_PEAK_CHAINS = [
    SegmentChain(GenomicSegment("chrI",100,15000,"+"),ID='feature1',signalValue=341.2,pValue=-123.2,qValue=-513.3,peak=50,_bedx_column_order=["signalValue","pValue","qValue","peak"],thickstart=100,thickend=100),
    SegmentChain(GenomicSegment("chrII",320,15000,"-"),ID='feature2',signalValue=2.1,pValue=-5123.2,qValue=0.0,peak=650,_bedx_column_order=["signalValue","pValue","qValue","peak"],thickstart=320,thickend=320),
]


