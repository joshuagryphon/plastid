#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.readers.wiggle`"""
import unittest
from plastid.util.services.mini2to3 import cStringIO
from nose.plugins.attrib import attr
from plastid.readers.wiggle import WiggleReader

#===============================================================================
# INDEX: test suites
#===============================================================================

@attr(test="unit")
class TestWiggleReader(unittest.TestCase):

    def _do(self,data_str,data_format,expected_results):
        """Execute tests on various formats

        Parameters
        ----------
        data_str : str
            wiggle file data

        data_format : str
            name of expected data format
        """
        reader = WiggleReader(cStringIO.StringIO(data_str))
        for n,(tup1, tup2) in enumerate(zip(expected_results,reader)):
            self.assertEquals(tup1,tup2)

        self.assertEquals(n,len(expected_results)-1,"Not all results processed: %s vs %s" % (n,len(expected_results)))
        
        self.assertEquals(reader.data_format,data_format)

    def test_read_multispan_multistep_varstep(self):
        self._do(_MULTISPAN_VARSTEP,"variableStep",_MULTISPAN_TUPLES)

    def test_read_multispan_multistep_fixedstep(self):
        self._do(_MULTISPAN_FIXEDSTEP,"fixedStep",_MULTISPAN_TUPLES)

    def test_read_multispan_multistep_bedgraph(self):
        self._do(_MULTISPAN_BEDGRAPH,"bedGraph",_MULTISPAN_TUPLES)


#===============================================================================
# INDEX: test data
#===============================================================================

# expected results from each type of reader
_MULTISPAN_TUPLES = [
                 ("chrA",5,10,5),
                 ("chrA",10,15,32),
                 ("chrA",25,30,40),
                 ("chrA",2000,2005,50),
                 ("chrA",2005,2010,50.5),
                 ("chrA",2010,2015,50.2),
                 ("chrB",0,20,1),
                 ("chrB",20,40,1),
                 ("chrB",1000,1020,30),
                 ("chrB",1020,1040,60),
                 ("chrB",1040,1060,1),
                 ("chrB",1060,1062,5),
                 ("chrB",1070,1072,10),
                 ("chrB",1080,1082,20),
                 ]

_MULTISPAN_BEDGRAPH="""track type=bedGraph
chrA    5   10  5
chrA    10  15  32
chrA    25  30  40
chrA    2000    2005    50
chrA    2005    2010    50.5
chrA    2010    2015    50.2
chrB    0   20  1
chrB    20  40  1
chrB    1000    1020    30
chrB    1020    1040    60
chrB    1040    1060    1
chrB    1060    1062    5
chrB    1070    1072    10
chrB    1080    1082    20
""".replace("    ","\t")

_MULTISPAN_VARSTEP="""track type=wiggle_0
variableStep chrom=chrA span=5
6   5
11  32
26  40
2001    50
2006    50.5
2011    50.2
variableStep chrom=chrB span=20
1   1
21  1
1001    30
1021    60
1041    1
variableStep chrom=chrB span=2
1061    5
1071    10
1081    20
""".replace("    ","\t")

_MULTISPAN_FIXEDSTEP="""track type=wiggle_0
fixedStep chrom=chrA start=6 span=5 step=5
5
32
fixedStep chrom=chrA start=26 span=5 step=5
40
fixedStep chrom=chrA start=2001 span=5 step=5
50
50.5
50.2
fixedStep chrom=chrB span=20 step=20
1
1
fixedStep chrom=chrB start=1001 span=20 step=20
30
60
1
fixedStep chrom=chrB start=1061 span=2 step=10
5
10
20
"""
