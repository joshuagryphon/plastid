#!/usr/bin/env python
"""Test suite for py:mod:`plastid.util.services.misc`"""
import numpy
import unittest
from nose.plugins.attrib import attr
from plastid.util.services.misc import guess_formatter, number

@attr(test="unit")
class TestMisc(unittest.TestCase):
    def setUp(self):
        self.tests = [("nan",numpy.nan),
                    ("Nan",numpy.nan),
                    ("NaN",numpy.nan),
                    ("None",numpy.nan),
                    ("none",numpy.nan),
                    ("Inf",numpy.inf),
                    ("inf",numpy.inf),
                    ("-Inf",-numpy.inf),
                    ("-inf",-numpy.inf),
                    ("5",5),
                    ("5.0",5),
                    ("-5",-5),
                    ("-5.0",-5),
                    ("0",0),
                    ("5.1",5.1),
                    ("-5.1",-5.1),
                    ("a5","a5"),
                    ("5a","5a"),
                    ("-5.1a","-5.1a"),
                    ("some_string","some_string"),
                    ("True",True),
                    ("False",False),
                    ("true",True),
                    ("false",False),
                    ("True5","True5"),
                    ("5True","5True"),
                    ("1e10",1e10),
                    ("1e-10",1e-10),
                    ("-1e10",-1e10),
                    ("-1e-10",-1e-10),
                   ]

    def test_guess_formatter(self):
        """Test parsing of strings to floats, ints, nans, infs, bools, and strings"""
        for source,dest in self.tests:
            if isinstance(dest,float):
                if not numpy.isnan(dest):
                    self.assertEquals(guess_formatter(source),dest)
                else:
                    self.assertTrue(numpy.isnan(guess_formatter(source)))
            else: #if isinstance(dest,int):
                self.assertEquals(guess_formatter(source),dest)

    def test_number(self):
        """Test numerical parsing of strings to floats, ints, nans, and infs"""
        for source,dest in self.tests:
            if isinstance(dest,bool):
                pass
            elif isinstance(dest,float):
                if not numpy.isnan(dest):
                    self.assertEquals(number(source),dest)
                else:
                    self.assertTrue(numpy.isnan(number(source)))
            elif isinstance(dest,int):
                self.assertEquals(number(source),dest)
            else:
                self.assertRaises(ValueError,number,source)
