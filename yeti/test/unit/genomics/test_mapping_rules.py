#!/usr/bin/env python
import numpy
import mock
import warnings
from nose.tools import assert_true, assert_equal
from yeti.genomics.genome_array import FivePrimeMapFactory,\
                                       ThreePrimeMapFactory,\
                                       CenterMapFactory,\
                                       VariableFivePrimeMapFactory
from yeti.genomics.roitools import GenomicSegment


class MockRead(object):
    """Mock class to substitute pysam.AlignedSegment"""
    def __init__(self,positions,strand):
        self.positions = positions
        self.is_reverse = True if strand == "+" else False


class TestBAMMappingRules(object):

    @classmethod
    def setUpClass(cls):
        min_ = 25
        max_ = 40
        cls.strands  = ("+","-")
        cls.segs = { X : GenomicSegment("mock",0,2000,X) for X in cls.strands }
        cls.reads = { Y : [MockRead(range(X),Y) for X in range(min_,max_)] for Y in ("+","-") }
        cls.expected = {}
        for mapping in ("fiveprime","threeprime","center"):
            for param in (0,10):
                for strand in cls.strands:
                    cls.expected[(mapping,param,strand)] = numpy.zeros(2000)

        cls.expected[("fiveprime",0, "+")][0]  = max_ - min_
        cls.expected[("fiveprime",10,"+")][10] = max_ - min_
        cls.expected[("fiveprime",0, "-")][min_-1:max_-1]   = 1
        cls.expected[("fiveprime",10,"-")][min_-11:max_-11] = 1

        cls.expected[("threeprime",0, "-")][0]  = max_ - min_
        cls.expected[("threeprime",10,"-")][10] = max_ - min_
        cls.expected[("threeprime",0, "+")][min_-1:max_-1]   = 1
        cls.expected[("threeprime",10,"+")][min_-11:max_-11] = 1

        for my_len in range(min_,max_):
            cls.expected[("center",0,"+")][:my_len] += 1.0/my_len
            cls.expected[("center",0,"-")][:my_len] += 1.0/my_len
            cls.expected[("center",10,"+")][10:my_len-10] += 1.0/(my_len-2*10)
            cls.expected[("center",10,"-")][10:my_len-10] += 1.0/(my_len-2*10)

        cls.map_factories = {
            "fiveprime"          : FivePrimeMapFactory,
            "threeprime"         : ThreePrimeMapFactory,
            "fiveprime_variable" : VariableFivePrimeMapFactory,
            "center"             : CenterMapFactory
        }

    def check_mapping(self,map_name,map_param,strand):
        fn = self.map_factories[map_name](map_param)
        expected = self.expected[(map_name,map_param,strand)]
        reads_out, count_array = self.map_factories[map_name](map_param)(self.reads[strand],self.segs[strand])
        assert_equal(reads_out,self.reads[strand])
        assert_true((count_array == expected).all())

    def test_fiveprime_threeprime_center(self):
        for mapping in ("fiveprime","threeprime","center"):
            for offset in 0,10:
                for strand in self.strands:
                    yield self.check_mapping, mapping, offset, strand

    def test_fiveprime_variable(self):
        offset_dicts = {}
        offset_dicts[("default_only","+")] = ({ "default" : 0 },
                                               self.expected[("fiveprime",0,"+")],
                                             )
        for (name, strand), (dict_, expected) in offset_dicts.items():
            fn = VariableFivePrimeMapFactory(dict_)
            reads_out, count_array = fn(self.reads[strand],self.segs[strand])
            assert_true((count_array==expected).all())

    def test_unmappable_raises_warnings(self):
        assert False

    def test_unmappable_not_mapped_in_vector_or_returned(self):
        assert False
