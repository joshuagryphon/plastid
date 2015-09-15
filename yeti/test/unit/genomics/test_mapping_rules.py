#!/usr/bin/env python
import numpy
import pysam
import warnings
from nose.tools import assert_true, assert_equal, assert_greater_equal
from plastid.genomics.map_factories import FivePrimeMapFactory,\
                                       ThreePrimeMapFactory,\
                                       CenterMapFactory,\
                                       VariableFivePrimeMapFactory
from plastid.genomics.roitools import GenomicSegment


class TestBAM_MappingRules(object):

    @classmethod
    def setUpClass(cls):
        min_ = 25
        max_ = 40
        cls.strands  = ("+","-")
        cls.segs = { X : GenomicSegment("mock",0,2000,X) for X in cls.strands }
        
        cls.reads = { Y : [cls.make_alignment(0,X,Y) for X in range(min_,max_)] for Y in ("+","-") }
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

    @staticmethod
    def make_alignment(start_pos,end_pos,strand):
        read = pysam.AlignedSegment()
        read.reference_start = start_pos
        my_length = end_pos - start_pos
        read.query_sequence = "N"*my_length
        read.reference_id = 0
        read.is_reverse = True if strand == "-" else False
        read.cigarstring = "%sM" % my_length
        return read

    def check_mapping(self,map_name,map_param,strand):
        fn = self.map_factories[map_name](map_param)
        expected = self.expected[(map_name,map_param,strand)]
        reads_out, count_array = self.map_factories[map_name](map_param)(self.reads[strand],self.segs[strand])
        msg1 = "failed to return reads with %s mapping with param %s on strand '%s'." % (map_name, map_param, strand)
        msg2 = "failed %s mapping vector with param %s on strand '%s'." % (map_name, map_param, strand)
        assert_equal(reads_out,self.reads[strand],msg1)
        assert_true((count_array == expected).all(),msg2)

    def test_fiveprime_threeprime_center(self):
        for mapping in ("fiveprime","threeprime","center"):
            for offset in 0,10:
                for strand in self.strands:
                    yield self.check_mapping, mapping, offset, strand

    def test_fiveprime_variable(self):
        offset_dicts = {}
        fancy = { X : X//2 for X in range (25,40) }
        fancy_reads = { X : numpy.zeros(2000) for X in self.strands }
        for read in self.reads["+"]:
            fancy_reads["+"][read.positions[fancy[len(read.positions)]]] += 1
        for read in self.reads["-"]:
            fancy_reads["-"][read.positions[-fancy[len(read.positions)]-1]] += 1

        fancy_plus  = numpy.zeros(2000)
        fancy_minus = numpy.zeros(2000)


        offset_dicts[("default_only","+")] = ({ "default" : 0 },
                                               self.expected[("fiveprime",0,"+")])
        offset_dicts[("default_only","-")] = ({ "default" : 0 },
                                               self.expected[("fiveprime",0,"-")])
        offset_dicts[("fancy","+")       ] = (fancy,fancy_reads["+"])
        offset_dicts[("fancy","-")       ] = (fancy,fancy_reads["-"])
                                             
        for (name, strand), (dict_, expected) in offset_dicts.items():
            fn = VariableFivePrimeMapFactory(dict_)
            reads_out, count_array = fn(self.reads[strand],self.segs[strand])
            msg = "Failed fiveprime variable mapping for test %s(%s)" % (name,strand)
            assert_true((count_array==expected).all(),msg)

    def check_unmappable_raises_warnings(self,test_name,map_param,strand):
        map_factory = self.map_factories[test_name]
        msg = "Map function '%s' failed to raise warning for unmappable reads on strand '%s'" % (test_name,strand)
        with warnings.catch_warnings(record=True) as warns:
            warnings.simplefilter("always")
            fn = map_factory(map_param)
            reads_out, count_array = fn(self.reads[strand],self.segs[strand])

        assert_greater_equal(len(warns),0,msg)

    def test_unmappable_raises_warnings(self):
        func_params = { "fiveprime"  : 30,
                        "threeprime" : 30,
                        "center"     : 15,
                        "fiveprime_variable" : { 25 : 10, "default" : 28 }
                       }
        for test_name, map_param in func_params.items():
            for strand in self.strands:
                yield self.check_unmappable_raises_warnings, test_name, map_param, strand 

    def check_unmappable_not_mapped_in_vector_or_returned(self,test_name,map_param,expected_num_reads,strand):
        msg1 = "Map function '%s' failed to filter out short reads on strand '%s' from returned reads" % (test_name,strand)
        msg2 = "Map function '%s' failed to filter out short reads on strand '%s' from count vector" % (test_name,strand)
        fn = self.map_factories[test_name](map_param)
        reads_out, count_array = fn(self.reads[strand],self.segs[strand])
        assert_equal(len(reads_out),expected_num_reads,msg1)
        assert_equal(count_array.sum(),expected_num_reads,msg2)

    def test_unmappable_not_mapped_in_vector_or_returned(self):
        func_params = { "fiveprime"  : 30,
                        "threeprime" : 30,
                        "center"     : 15,
                        "fiveprime_variable" : { 25 : 10, "default" : 28 }
                       }
        expected_num = { "fiveprime"  : len([X for X in self.reads["+"] if len(X.positions) > func_params["fiveprime"]]),
                         "threeprime" : len([X for X in self.reads["+"] if len(X.positions) > func_params["threeprime"]]),
                         "center"     : len([X for X in self.reads["+"] if len(X.positions) > 2*func_params["center"]]),
                         "fiveprime_variable" : len([X for X in self.reads["+"] if len(X.positions) > 28 or len(X.positions) == 25]),
                       }
        for test_name, map_param in func_params.items():
            expected_num_reads = expected_num[test_name]
            for strand in self.strands:
                yield self.check_unmappable_not_mapped_in_vector_or_returned, test_name, map_param, expected_num_reads, strand 

