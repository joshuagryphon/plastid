#!/usr/bin/env python
"""
"""
import unittest
import copy
import numpy
import pandas as pd
from random import shuffle
from pkg_resources import resource_filename, cleanup_resources
from nose.plugins.attrib import attr
from plastid.genomics.roitools import GenomicSegment


# slight hack to keep imported method from being run as a test
# can't use unittest.skip, or else no tests will never be run!
from plastid.bin.test_table_equality import test_dataframe_equality as checkeq
checkeq.__name__   = "checkeq"
checkeq.__module__ = "checkeq"

# components we will use in equality tests
size = 5000


@attr(test="unit")
class TestTestDataframeEquality(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.cols = { "intA"   : numpy.random.randint(0,high=2**16,size=size),
         "intB"   : numpy.random.randint(-10,high=20,size=size),
         "idxA"   : numpy.arange(size),
         "chrA"   : numpy.array([chr(65+(X%(91-65))) for X in range(size)]),
         "strA"   : numpy.array([str(GenomicSegment("chrA",X,X+500,"+")) for X in range(size)]),
         "strB"   : numpy.array([str(GenomicSegment("chrB",X/2,X/2+500,"-")) for X in range(size)]),
         "floatA" : 10*numpy.random.randn(size) + 500,
         "floatB" : (10**-5)*numpy.random.random(size),
         "objA"   : numpy.tile(None,5000),
         "objB"   : numpy.array([GenomicSegment("chrC",X,X+Y,"+") for X,Y in zip(range(size),numpy.random.randint(2,high=1000,size=size))]),
       }

    def test_dataframe_quality_when_identical(self):
        df1 = pd.DataFrame(self.cols)
        self.assertTrue(checkeq(df1,df1))

    def test_dataframe_equality_no_sort(self):
        df1 = pd.DataFrame(self.cols)
        df2 = pd.DataFrame(copy.deepcopy(self.cols))
        self.assertTrue(checkeq(df1,df2))
        self.assertTrue(checkeq(df2,df1))
    
    def test_dataframe_equality_within_tol(self):
        tol = 10**-8
        noiseA = tol/10**2 * numpy.random.randn(size)
        noiseB = tol/10**2 * numpy.random.randn(size)
        df1 = pd.DataFrame(self.cols)
        df2 = copy.deepcopy(df1)
        df2["floatA"] += noiseA
        df2["floatB"] -= noiseB
        self.assertTrue(checkeq(df1,df2,tol=tol))
        self.assertTrue(checkeq(df2,df1,tol=tol))

    def test_dataframe_inequality_above_tol(self):
        tol = 10**-8
        noiseA = tol*10**2 * numpy.random.randn(size)
        noiseB = tol*10**2 * numpy.random.randn(size)
        df1 = pd.DataFrame(self.cols)
        df2 = copy.deepcopy(df1)
        df2["floatA"] += noiseA
        df2["floatB"] -= noiseB
        self.assertFalse(checkeq(df1,df2,tol=tol))
        self.assertFalse(checkeq(df2,df1,tol=tol))

    def test_dataframe_inequality_wrong_columns(self):
        df1 = pd.DataFrame(self.cols)
        df2 = pd.DataFrame({ K : copy.deepcopy(self.cols[K]) for K in sorted(self.cols.keys())[:-2] })
        self.assertFalse(checkeq(df1,df2))
        self.assertFalse(checkeq(df2,df1))
    
    def test_dataframe_inequality_wrong_rows(self):
        df1 = pd.DataFrame(self.cols)
        df2 = pd.DataFrame({K : self.cols[K][:size-1000] for K in self.cols.keys()})
        self.assertFalse(checkeq(df1,df2))
        self.assertFalse(checkeq(df2,df1))
    
    def check_dataframe_equality_same(self,special_val):
        """Helper function for checking equality when various special values
        are in the same location in both dataframes

        Parameters
        ----------
        special_val : numpy.nan, numpy.inf, -numpy.inf, or None
            Value that is ignored by :py:meth:`plastid.bin.test_dataframe_equality` if it occurs in the same cells of both dataframes
        """
        idx = numpy.random.randint(0,high=size,size=500)
        tmpcols = copy.deepcopy(self.cols)
        tmpcols["floatA"][idx] = special_val
        df1 = pd.DataFrame(tmpcols)
        df2 = pd.DataFrame(copy.deepcopy(tmpcols))
        self.assertTrue(checkeq(df1,df2))
        self.assertTrue(checkeq(df2,df1))

    def test_dataframe_equality_same_nans(self):
        self.check_dataframe_equality_same(numpy.nan)

    def test_dataframe_equality_same_infs(self):
        self.check_dataframe_equality_same(numpy.inf)
        self.check_dataframe_equality_same(-numpy.inf)

    def check_dataframe_inequality_different(self,special_val):
        """Helper function for checking inequality when various special values
        appear in the different locations in both dataframes

        Parameters
        ----------
        special_val : numpy.nan, numpy.inf, -numpy.inf, or None
            Value that is ignored by :py:meth:`plastid.bin.test_dataframe_equality` if it occurs in the same cells of both dataframes
        """
        idxA = numpy.random.randint(0,high=size,size=500)
        idxB = numpy.random.randint(0,high=size,size=500)

        colsA = copy.deepcopy(self.cols)
        colsA["floatA"][idxA] = special_val

        colsB = copy.deepcopy(self.cols)
        colsB["floatA"][idxB] = special_val

        df1 = pd.DataFrame(colsA)
        df2 = pd.DataFrame(colsB)
        self.assertFalse(checkeq(df1,df2))
        self.assertFalse(checkeq(df2,df1))

    def test_dataframe_inequality_different_nans(self):
        self.check_dataframe_inequality_different(numpy.nan)

    def test_dataframe_inequailty_different_infs(self):
        self.check_dataframe_inequality_different(numpy.inf)
        self.check_dataframe_inequality_different(-numpy.inf)

    def test_dataframe_equality_with_sort_numeric(self):
        shuffidx = numpy.arange(size)
        shuffle(shuffidx)
        df1 = pd.DataFrame(self.cols)
        df2 = pd.DataFrame({ K : self.cols[K][shuffidx] for K in self.cols.keys()})
        self.assertTrue(checkeq(df1,df2,sort_columns=["strA"]))
        self.assertTrue(checkeq(df2,df1,sort_columns=["strA"]))
    
    def test_dataframe_equality_with_sort_str(self):
        shuffidx = numpy.arange(size)
        shuffle(shuffidx)
        df1 = pd.DataFrame(self.cols)
        df2 = pd.DataFrame({ K : self.cols[K][shuffidx] for K in self.cols.keys()})
        self.assertTrue(checkeq(df1,df2,sort_columns=["idxA"]))
        self.assertTrue(checkeq(df2,df1,sort_columns=["idxA"]))

    def test_dataframe_equality_without_sort(self):
        shuffidx = numpy.arange(size)
        shuffle(shuffidx)
        df1 = pd.DataFrame(self.cols)
        df2 = pd.DataFrame({ K : self.cols[K][shuffidx] for K in self.cols.keys()})
        self.assertFalse(checkeq(df1,df2))
        self.assertFalse(checkeq(df2,df1))
 
    def test_dataframe_equality_with_multi_sort(self):
        shuffidx = numpy.arange(size)
        shuffle(shuffidx)
        df1 = pd.DataFrame(self.cols)
        df2 = pd.DataFrame({ K : self.cols[K][shuffidx] for K in self.cols.keys()})
        self.assertTrue(checkeq(df1,df2,sort_columns=["strB","chrA"]))
        self.assertTrue(checkeq(df2,df1,sort_columns=["strB","chrA"]))
