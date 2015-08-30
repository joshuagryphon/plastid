#!/usr/bin/env python
"""Test suite for :py:mod:`yeti.readers.bed`

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
import pandas as pd

from nose.plugins.attrib import attr
from yeti.util.services.mini2to3 import cStringIO
from yeti.genomics.roitools import SegmentChain, GenomicSegment, Transcript
from yeti.readers.bed import BED_to_Transcripts, BED_to_SegmentChain, BED_Reader
from nose.tools import assert_equal

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
        cls.data[12] = _BED12_DATA
        bed_df = pd.read_table(cStringIO.StringIO(_BED12_DATA),header=None,sep="\t",index_col=None)
        extra_df = pd.read_table(cStringIO.StringIO(_EXTRA_COLS),header=0,sep="\t",index_col=None)
        cls.big_df = pd.concat([bed_df,extra_df],axis=1)

        for n in (3,4,5,6,8,9):
            #cls.data[n] = cls.get_bed_subset(cls.header,cls.data[12],n)
            cls.data[n] = cls.get_bed_subset2(cls.header,n,0)

    @classmethod
    def get_bed_subset2(cls,header,bed_cols,extra_cols=0):
        buf = cStringIO.StringIO()
        columns = cls.big_df.columns[list(range(bed_cols)) + list(range(bed_cols,bed_cols+extra_cols))]
        cls.big_df.to_csv(buf,columns=columns,sep="\t",index=False,header=False,float_format="%.8f")
        return buf.getvalue()
        
    @staticmethod
    def get_bed_subset(header,data,columns):
        """Select a subset of columns from BED data block
        
        Parameters
        ----------
        header : str
            Multi-line header data for BED file
            
        data : str
            Multi-line BED12 data
        
        columns : int
            Number of columns to include
        
        Returns
        -------
        str : block of BED data, including header
        """
        ltmp = []
        for line in data.split("\n"):
            ltmp.append("%s\n" % "\t".join(line.split("\t")[:columns]))
        
        return header+"".join(ltmp)
   
    @staticmethod
    def check_equal(found,expected,msg=None):
        if msg is not None:
            assert_equal(found,expected,msg)
        else:
            assert_equal(found,expected)

    def test_bed_to_various(self):
        """Helper function for BED import to SegmentChain or Transcripts
        
        Parameters
        ----------
        known_set : list<SegmentChain> or list<Transcript>
            Expected results of parsing BED blocks
        
        reader_fn : Function or class
            Parser for blocks of BED formatted data
        """
        tx_reader = functools.partial(BED_Reader,return_type=Transcript)
        tests = [(BED_to_SegmentChain,_TEST_IVCOLLECTIONS,"tosegmentchain,segmentchain"),
                 (BED_to_Transcripts,_TEST_TRANSCRIPTS,"totranscripts_transcript"),
                 (BED_Reader,_TEST_IVCOLLECTIONS,"reader_segmentchain"),
                 (tx_reader,_TEST_TRANSCRIPTS,"reader_transcript"),
                ]
        for reader_fn, known_set, name in tests: #BED_to_SegmentChain, BED_to_Transcripts:
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

    def test_ivcollection_thick_start_end(self):
        """Checks equality of thickstart and thickend attributes for SegmentChain objects"""
        for n,data_str in sorted(self.data.items()):
            for c, (test_ivc,known_ivc) in enumerate(zip(BED_to_SegmentChain(cStringIO.StringIO(data_str)),
                                                       _TEST_IVCOLLECTIONS)):
                if n >= 8:
                    err_msg = "Failed thickstart/end equality on %s-column BED input: %s,%s" % (n,known_ivc.attr,test_ivc.attr)
                    if known_ivc.attr.get("thickstart",None) is not None:
                        yield self.check_equal, known_ivc.attr["thickstart"],test_ivc.attr["thickstart"],err_msg
                    if known_ivc.attr.get("thickend",None) is not None:
                        yield self.check_equal, known_ivc.attr.get("thickend"),test_ivc.attr["thickend"],err_msg
        
            yield self.check_equal, c,20-1,"Not all intervals loaded! Expected %s, found %s." % (20-1,c)


    def test_transcript_cds_start_end_helper(self):
        """Checks equality of endpoints of coding regions for Transcript objects"""
        for n,data_str in sorted(self.data.items()):
            for c, (test_ivc,known_ivc) in enumerate(zip(BED_to_Transcripts(cStringIO.StringIO(data_str)),
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
        

#===============================================================================
# INDEX: test data
#===============================================================================

# test dataset, constructed manually to include various edge cases
_TEST_IVCOLLECTIONS = [
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
_TEST_TRANSCRIPTS = [Transcript(*X._segments,**X.attr) for X in _TEST_IVCOLLECTIONS]


_BED_HEADER = """browser position chrA:100-1100
track name=test_data description='my test data'
"""

# same data, as BED12 block
_BED12_DATA = """chrA    100    1100    IVC1p    0.0    +    -1    -1    0,0,0    1    1000,    0,
chrA    100    1100    IVC1m    0.0    -    -1    -1    0,0,0    1    1000,    0,
chrA    100    2600    IVC2p    0.0    +    -1    -1    0,0,0    2    1000,500,    0,2000,
chrA    100    2600    IVC2m    0.0    -    -1    -1    0,0,0    2    1000,500,    0,2000,
chrA    100    2600    IVC3p    500.0    +    -1    -1    0,0,0    2    1000,500,    0,2000,
chrA    100    2600    IVC3m    500.0    -    -1    -1    0,0,0    2    1000,500,    0,2000,
chrA    100    2700    IVC4p    500.0    +    -1    -1    0,0,0    3    1000,500,95,    0,2000,2505,
chrA    100    2600    IVC4m    500.0    -    -1    -1    0,0,0    2    1000,500,    0,2000,
chrA    100    2700    IVC5p    500.0    +    -1    -1    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC5m    500.0    -    -1    -1    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC6p    500.0    +    -1    -1    0,122,223    3    1000,500,95,    0,2000,2505,
chrA    100    2700    IVC6m    500.0    -    -1    -1    0,122,223    3    1000,500,95,    0,2000,2505,
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
2    30.1235    DEF    gene_id "gene_2"; transcript_id "transcript_2";
3    15123.2    ghi    gene_id "gene_3"; transcript_id "transcript_3";
4    2.0    alongword    gene_id "gene_4"; transcript_id "transcript_4";
5    -3.1234    a sentence with spaces    gene_id "gene_5"; transcript_id "transcript_5";
6    -20.5    some notes with "quotes"    gene_id "gene_6"; transcript_id "transcript_6";
7    -1e10    1    gene_id "gene_7"; transcript_id "transcript_7";
8    2e5    2    gene_id "gene_8"; transcript_id "transcript_8";
9    2.3e6    3.0    gene_id "gene_9"; transcript_id "transcript_9";
10    0.000003    string1    gene_id "gene_10"; transcript_id "transcript_10";
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
