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
from yeti.util.services.mini2to3 import cStringIO
import functools
from nose.plugins.attrib import attr
from yeti.genomics.roitools import SegmentChain, GenomicSegment, Transcript
from yeti.readers.bed import BED_to_Transcripts, BED_to_SegmentChain, BED_Reader

#===============================================================================
# INDEX: test suites
#===============================================================================
@attr(test="unit")
class TestBED(unittest.TestCase):
    """Test case for BED input/output"""
    
    @classmethod
    def setUpClass(cls):
        cls.header = _BED_HEADER
        cls.data = {}
        cls.data[12] = _BED12_DATA
        for n in (3,4,5,6,8,9):
            cls.data[n] = cls.get_bed_subset(cls.header,cls.data[12],n)
        
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
    
    def bed_to_various_helper(self,known_set,reader_fn):
        """Helper function for BED import to SegmentChain or Transcripts
        
        Parameters
        ----------
        known_set : list<SegmentChain> or list<Transcript>
            Expected results of parsing BED blocks
        
        reader_fn : Function or class
            Parser for blocks of BED formatted data
        """
        for n,data_str in sorted(self.data.items()):
            c = 0
            for (test_ivc,known_ivc) in zip(reader_fn(cStringIO.StringIO(data_str)),
                                                       known_set):
                # columns: chrom, start, end
                if n >= 3:
                    # no strand info, so we need to test iv.start, iv.end, iv.chrom
                    err_msg = "Failed endpoint equality on %s-column BED input: %s,%s" % (n,known_ivc,test_ivc)
                    self.assertEqual(known_ivc.spanning_segment.start,test_ivc.spanning_segment.start,err_msg)
                    self.assertEqual(known_ivc.spanning_segment.end,test_ivc.spanning_segment.end,err_msg)
                    self.assertEqual(known_ivc.spanning_segment.chrom,test_ivc.spanning_segment.chrom,err_msg)
                # column: name
                if n >= 4:
                    err_msg = "Failed name equality on %s-column BED input: %s,%s" % (n,known_ivc.attr,test_ivc.attr)
                    self.assertEqual(known_ivc.attr["ID"],test_ivc.attr["ID"],err_msg)
                # column: score
                if n >= 5:
                    err_msg = "Failed score equality on %s-column BED input: %s,%s" % (n,known_ivc.attr,test_ivc.attr)
                    self.assertEqual(known_ivc.attr.get("score",0),test_ivc.attr["score"],err_msg)
                # column : strand
                if n >= 6:
                    err_msg = "Failed strand equality on %s-column BED input: %s,%s" % (n,known_ivc,test_ivc)
                    self.assertEqual(known_ivc.spanning_segment.strand,test_ivc.spanning_segment.strand)
                # column: color
                if n >= 9:
                    err_msg = "Failed color equality on %s-column BED input: %s,%s" % (n,known_ivc.attr,test_ivc.attr)
                    self.assertEqual(known_ivc.attr.get("color","#000000"),test_ivc.attr["color"],err_msg)
                # columns: exon/block info
                if n == 12:
                    err_msg = "Failed block equality on %s-column BED input: %s,%s" % (n,known_ivc,test_ivc)
                    for iv1, iv2 in zip(known_ivc,test_ivc):
                        self.assertEqual(iv1,iv2,err_msg)
                    err_msg = "Failed position set on %s-column BED input: %s,%s" % (n,known_ivc,test_ivc)
                    self.assertEqual(known_ivc.get_position_set(),test_ivc.get_position_set(),err_msg)
                
                c += 1
            
            self.assertEqual(c,len(known_set),"Not all intervals loaded! Expected %s, found %s." % (len(known_set),c))

    def ivcollection_thick_start_end_helper(self):
        """Checks equality of thickstart and thickend attributes for SegmentChain objects"""
        for n,data_str in sorted(self.data.items()):
            for c, (test_ivc,known_ivc) in enumerate(zip(BED_to_SegmentChain(cStringIO.StringIO(data_str)),
                                                       _TEST_IVCOLLECTIONS)):
                if n >= 8:
                    err_msg = "Failed thickstart/end equality on %s-column BED input: %s,%s" % (n,known_ivc.attr,test_ivc.attr)
                    if known_ivc.attr.get("thickstart",None) is not None:
                        self.assertEqual(known_ivc.attr["thickstart"],test_ivc.attr["thickstart"],err_msg)
                    if known_ivc.attr.get("thickend",None) is not None:
                        self.assertEqual(known_ivc.attr.get("thickend"),test_ivc.attr["thickend"],err_msg)
        
            self.assertEqual(c,20-1,"Not all intervals loaded! Expected %s, found %s." % (20-1,c))


    def transcript_cds_start_end_helper(self):
        """Checks equality of endpoints of coding regions for Transcript objects"""
        for n,data_str in sorted(self.data.items()):
            for c, (test_ivc,known_ivc) in enumerate(zip(BED_to_Transcripts(cStringIO.StringIO(data_str)),
                                                       _TEST_TRANSCRIPTS)):
                if n >= 8:
                    err_msg = "Failed thickstart/end equality on %s-column BED input: %s,%s" % (n,known_ivc.attr,test_ivc.attr)
                    if known_ivc.attr.get("cds_genome_start",None) is not None:
                        self.assertEqual(known_ivc.attr["cds_start"],test_ivc.attr["cds_start"],err_msg)
                        self.assertEqual(known_ivc.attr["cds_genome_start"],test_ivc.attr["cds_genome_start"],err_msg)
                        self.assertEqual(known_ivc.cds_genome_start,test_ivc.cds_genome_start,err_msg)
                        self.assertEqual(known_ivc.cds_start,test_ivc.cds_start,err_msg)
                    if known_ivc.attr.get("cds_genome_end",None) is not None:
                        self.assertEqual(known_ivc.attr["cds_end"],test_ivc.attr["cds_end"],err_msg)
                        self.assertEqual(known_ivc.attr["cds_genome_end"],test_ivc.attr["cds_genome_end"],err_msg)
                        self.assertEqual(known_ivc.cds_genome_end,test_ivc.cds_genome_end,err_msg)
                        self.assertEqual(known_ivc.cds_end,test_ivc.cds_end,err_msg)
            
            self.assertEqual(c,20-1,"Not all intervals loaded! Expected %s, found %s." % (20-1,c))
        
    def test_bed_to_ivcollection(self):
        self.bed_to_various_helper(reader_fn=BED_to_SegmentChain,known_set=_TEST_IVCOLLECTIONS)
        self.ivcollection_thick_start_end_helper()

    def test_bed_to_transcripts(self):
        self.bed_to_various_helper(reader_fn=BED_to_Transcripts,known_set=_TEST_TRANSCRIPTS)
        self.transcript_cds_start_end_helper()
    
    def test_bed_reader_to_ivcollections(self):
        self.bed_to_various_helper(reader_fn=BED_Reader,known_set=_TEST_IVCOLLECTIONS)
        self.ivcollection_thick_start_end_helper()

    def test_bed_reader_to_transcripts(self):
        partial = functools.partial(BED_Reader,return_type=Transcript)
        self.bed_to_various_helper(reader_fn=partial,known_set=_TEST_TRANSCRIPTS)
        self.transcript_cds_start_end_helper()

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
