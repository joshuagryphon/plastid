#!/usr/bin/env python
"""Tests for data structures defined in :py:mod:`plastid.readers.common`
"""
import unittest
import copy
from nose.plugins.attrib import attr
from plastid.readers.common import add_three_for_stop_codon, \
                                                get_identical_attributes
from plastid.genomics.roitools import GenomicSegment,\
                                         SegmentChain, \
                                         Transcript

@attr(test="unit")
class TestAddThreeForStopCodon(unittest.TestCase):
    """Tests :py:func:`add_three_for_stop_codon`"""
    
    def test_noncoding_plus_strand(self):
        iv1p = GenomicSegment("chrA",100,190,"+")
        iv2p = GenomicSegment("chrA",200,203,"+")
        
        tx = add_three_for_stop_codon(Transcript(iv1p,iv2p))
        self.assertTrue(tx.cds_genome_start is None)
        self.assertTrue(tx.cds_start is None)
        self.assertTrue(tx.cds_genome_end is None)
        self.assertTrue(tx.cds_end is None)

    def test_noncoding_minus_strand(self):
        iv1m = GenomicSegment("chrA",100,190,"-")
        iv2m = GenomicSegment("chrA",200,203,"-")
        
        tx = add_three_for_stop_codon(Transcript(iv1m,iv2m))
        self.assertTrue(tx.cds_genome_start is None)
        self.assertTrue(tx.cds_start is None)
        self.assertTrue(tx.cds_genome_end is None)
        self.assertTrue(tx.cds_end is None)
                
    def test_no_autogrow_plus_strand(self):
        iv1p = GenomicSegment("chrA",100,190,"+")
        iv2p = GenomicSegment("chrA",200,203,"+")
        iv3p = GenomicSegment("chrA",200,201,"+")
        iv4p = GenomicSegment("chrA",204,206,"+")
        iv5p = GenomicSegment("chrA",220,250,"+")
        
        tests = [
            # new end is at end of only exon
            (Transcript(iv1p,cds_genome_start=100,cds_genome_end=187),190),
            
            # old end at end of first exon; new end is at end of final exon
            (Transcript(iv1p,iv2p,cds_genome_start=100,cds_genome_end=190),203),

            # old end at end of first exon; new end is at end of second of three exons
            (Transcript(iv1p,iv2p,iv3p,cds_genome_start=100,cds_genome_end=190),203),

            # old end at end of first exon; new end is at end of second of three exons
            # but end is listed in transcript-relative coordinates
            (Transcript(iv1p,iv2p,iv3p,cds_genome_start=100,cds_genome_end=200),203),

            
            # old end at end of first exon; new end at end of third/final exon
            (Transcript(iv1p,iv3p,iv4p,cds_genome_start=100,cds_genome_end=190),206),

            # old end at end of first exon; new end at end of third of four exons,
            # but position not actually in transcript
            #
            # stop codon split across exons 2 and 3
            (Transcript(iv1p,iv3p,iv4p,iv5p,cds_genome_start=100,cds_genome_end=190),206),

            
            # old end inside first exon; new end in third exon            
            (Transcript(iv1p,iv3p,iv4p,cds_genome_start=100,cds_genome_end=189),205),
            
            # old end inside first exon; new end at end of second exon
            # but position not actually in transcript
            # stop codon split across exons 1 and 2
            (Transcript(iv1p,iv3p,iv4p,cds_genome_start=100,cds_genome_end=188),201),       
            ]
        
        for tx, expected_end in tests:
            new_tx = add_three_for_stop_codon(tx)
            self.assertEquals(new_tx.cds_genome_end,expected_end)
            self.assertEquals(new_tx.cds_end,tx.cds_end+3)
    
    def test_autogrow_plus_strand(self):
        iv1p = GenomicSegment("chrA",100,190,"+")
        iv2p = GenomicSegment("chrA",200,203,"+")
        iv3p = GenomicSegment("chrA",200,201,"+")
        iv4p = GenomicSegment("chrA",204,206,"+")
        
        tests = [
            # new end will be 193. tx will need to grow by 3. unspliced
            (Transcript(iv1p,cds_genome_start=100,cds_genome_end=190),193),
            
            # new end will be 206, spliced
            (Transcript(iv1p,iv2p,cds_genome_start=100,cds_genome_end=203),206),
            
            # new end will be 206, spliced
            (Transcript(iv1p,iv3p,cds_genome_start=100,cds_genome_end=201),204),

            # new end will be 209, spliced in stop codon
            (Transcript(iv1p,iv3p,iv4p,cds_genome_start=100,cds_genome_end=206),209),            
            ]

        for tx, expected_end in tests:
            new_tx = add_three_for_stop_codon(tx)
            self.assertEquals(new_tx.cds_genome_end,expected_end)
            self.assertEquals(new_tx.cds_end,tx.cds_end+3)

    def test_no_autogrow_minus_strand(self):
        iv1m = GenomicSegment("chrA",100,190,"-")
        iv2m = GenomicSegment("chrA",90,93,"-")
        iv3m = GenomicSegment("chrA",90,92,"-")
        iv4m = GenomicSegment("chrA",80,81,"-")
        
        tests = [
            # new start is at start of only exon
            (Transcript(iv1m,cds_genome_start=103,cds_genome_end=190),100),
            
            # old start at start of second exon; new start is start of first exon
            (Transcript(iv2m,iv1m,cds_genome_start=100,cds_genome_end=190),90),
            
            # old start at start of third exon; new start at start of first exon
            (Transcript(iv4m,iv3m,iv1m,cds_genome_start=100,cds_genome_end=190),80),
            
            # old in third exon; new start in second exon
            (Transcript(iv4m,iv3m,iv1m,cds_genome_start=101,cds_genome_end=190),90),
            
            # old end in third exon; new end at beginning of third exon
            # but position not actually in transcript
            (Transcript(iv4m,iv3m,iv1m,cds_genome_start=102,cds_genome_end=190),91),            
            ]
        
        for tx, expected_start in tests:
            new_tx = add_three_for_stop_codon(tx)
            self.assertEquals(new_tx.cds_genome_start,expected_start)
            #self.assertEquals(new_tx.cds_end,tx.cds_end + 3)
            
    def test_autogrow_minus_strand(self):
        iv1m = GenomicSegment("chrA",100,190,"-")
        iv2m = GenomicSegment("chrA",90,93,"-")
        iv3m = GenomicSegment("chrA",90,92,"-")
        iv4m = GenomicSegment("chrA",90,91,"-")
        
        tests = [
            (Transcript(iv1m,cds_genome_start=100,cds_genome_end=190),97),
            (Transcript(iv2m,iv1m,cds_genome_start=90,cds_genome_end=190),87),
            (Transcript(iv3m,iv1m,cds_genome_start=90,cds_genome_end=190),87),
            (Transcript(iv4m,iv1m,cds_genome_start=90,cds_genome_end=190),87),            
            ]

        for tx, expected_start in tests:
            new_tx = add_three_for_stop_codon(tx)
            self.assertEquals(new_tx.cds_genome_start,expected_start)
            self.assertEquals(new_tx.cds_start,0)


@attr(test="unit")
class TestGetIdenticalAttributes(unittest.TestCase):
    """Tests :py:func:`get_identical_attributes`"""
    
    @classmethod
    def setUpClass(cls):
        cls.ivs = [GenomicSegment("chrA",100,190,"+"),
                   GenomicSegment("chrA",200,203,"+"),
                   GenomicSegment("chrA",200,201,"+"),
                   GenomicSegment("chrA",204,206,"+"),
                  ]
        cls.common_attr = dict(common1="common",common2="also common",common3="still common")
        cls.attrs = [dict(common_diff_val="unique_f1",unique_f1_key="something"),
                     dict(common_diff_val="unique_f2",unique_f2_key="something",unique_f2f3="something else"),
                     dict(common_diff_val="unique_f3",unique_f3_key="something",unique_f2f3="something else",unique_f3f4="f3 only"),
                     dict(common_diff_val="unique_f4",unique_f4_key="something",unique_f3f4="f4 only"),
                    ]
        for x in cls.attrs:
            x.update(cls.common_attr)
        
    def test_on_SegmentChain_no_exclude(self):
        features = [SegmentChain(self.ivs[n],**self.attrs[n]) for n in range(len(self.ivs))]
        common_plus_type = copy.deepcopy(self.common_attr)
        common_plus_type["type"] = "exon"
        self.assertEqual(get_identical_attributes(features),common_plus_type)
    
    def test_on_SegmentChain_exclude(self):
        features = [SegmentChain(self.ivs[n],**self.attrs[n]) for n in range(len(self.ivs))]
        self.assertEqual(get_identical_attributes(features,exclude=["type"]),self.common_attr)
        
