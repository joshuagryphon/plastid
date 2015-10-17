#!/usr/bin/env python
"""Test suite for py:mod:`plastid.plotting.colors`"""
import unittest
import numpy
from plastid.plotting.colors import get_rgb255, get_str_from_rgb255
from nose.plugins.attrib import attr

@attr(test="unit")
class TestColors(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # valid colors as hex strings and tuples
        cls.valid_colors = [("#000000",(0,0,0)),
                            ("#FFFFFF",(255,255,255)),
                            ("#ffffff",(255,255,255)),
                            ("#FF0000",(255,0,0)),
                            ("#00FF00",(0,255,0)),
                            ("#0000FF",(0,0,255)),
                            ("#6A5ACD",(106,90,205)),
                            ("#ee6a50",(238,106,80)),
                            ("#bf3EFF",(191,62,255)),
                            ("#a2CD5a",(162,205,90)),
                           ]

        # data that should cause problems
        cls.error_colors = [
                            ("#ZZZZZZ",(-5,10,255)), # out of range
                            ("#-23a5*",(10,3000,0)), # nonsense/out of range
                            ("#AA",(0,2)),           # too short
                            ("#AAAAAAAA",(5,5,5,5)), # too long
                            ("#GG0011",(0.5,3,2)) # float
                           ] 

    def test_get_rgb255_match(self):
        # test correct conversions
        for my_str, my_tup in self.valid_colors:
            self.assertEqual(tuple(get_rgb255(my_str)),my_tup)

    def test_get_rgb255_notmatch(self):
        # assert not incorrect conversions
        randidx = numpy.random.randint(0,high=len(self.valid_colors),size=len(self.valid_colors))
        for i,(my_str,_) in zip(randidx,self.valid_colors):
            if self.valid_colors[i][0].upper() != my_str.upper():
                self.assertNotEqual(tuple(get_rgb255(my_str)),
                                    self.valid_colors[i][0])

    def test_get_rgb255_knownfail(self):
        # assert invalid colors cause problems
        for my_str, _ in self.error_colors:
            self.assertRaises(ValueError,get_rgb255,my_str)

    def test_get_str_from_rgb255_match(self):
        # test correct conversions
        for my_str, my_tup in self.valid_colors:
            self.assertEqual(my_str.upper(),get_str_from_rgb255(my_tup).upper())

    def test_get_str_from_rgb255_notmatch(self):
        # assert not incorrect conversions
        randidx = numpy.random.randint(0,high=len(self.valid_colors),size=len(self.valid_colors))
        for i,(_,my_tup) in zip(randidx,self.valid_colors):
            if self.valid_colors[i][1] != my_tup:
                self.assertNotEqual(get_str_from_rgb255(my_tup),
                                    self.valid_colors[i][1])

    def test_get_str_from_rgb255_knownfail(self):
        # assert invalid colors cause problems
        for _,my_tup in self.error_colors:
            self.assertRaises(ValueError,get_str_from_rgb255,my_tup)

