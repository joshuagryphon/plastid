#!/usr/bin/env python
"""Test cases for :py:mod:`plastid.util.io.binary`"""
from plastid.util.io.binary import BinaryParserFactory, find_null_bytes
import unittest
import re
from nose.plugins.attrib import attr
from io import BytesIO #from plastid.util.services.mini2to3 import cStringIO

@attr(test="unit")
class TestBinaryParserFactory(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.factories = { "color" : BinaryParserFactory("ColorFactory","BBB",["r","g","b"])
                        }
        
        cls.results = { "color" : [dict(zip(list("rgb"),X)) for X in [(255,0  ,0  ),
                                                                      (255,255,0  ),
                                                                      (255,255,255),
                                                                      (255,0  ,255),
                                                                      ]] 
                      }

        cls.rawdata = { "color" : b"\xFF\x00\x00"+\
                                  b"\xFF\xFF\x00"+\
                                  b"\xFF\xFF\xFF"+\
                                  b"\xFF\x00\xFF" 
                      }
        
        cls.size_dict = { "x"  : 0,
                          "c"  : 1,
                          "b"  : 1,
                          "B"  : 1,
                          "?"  : 1,
                          "h"  : 2,
                          "H"  : 2,
                          "i"  : 4,
                          "I"  : 4,
                          "l"  : 4,
                          "L"  : 4,
                          "q"  : 8,
                          "Q"  : 8,
                          "f"  : 4,
                          "d"  : 8,
                          "s"  : 1,
                          "p"  : 1,
                        }
        
        cls.fmt_pat = re.compile(r"(?P<size>[0-9]+)*(?P<code>\w)")

    def call_helper(self,key):
        fh = BytesIO(self.rawdata[key])
        for result in self.results[key]:
            self.assertEqual(self.factories[key](fh),result)
        
        fh.close
                
    def test_call_color(self):
        self.call_helper("color")

@attr(test="unit")
def test_find_null_bytes():
    test = b"12345\x00123456789\x001234567890\x00"
    assert list(find_null_bytes(test)) == [5,15,26]