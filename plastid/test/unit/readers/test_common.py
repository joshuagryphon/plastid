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
        
