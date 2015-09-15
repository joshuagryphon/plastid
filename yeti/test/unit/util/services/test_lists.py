#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.util.services.lists`"""
import unittest
import numpy
from plastid.util.services.lists import parse_list, flatten_nested_lists_to_generator, flatten_nested_lists_to_list
from nose.plugins.attrib import attr

@attr(test="unit")
class TestLists(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # list that will be entirely flattened
        cls.tree1 = ["a",
                     "b",
                     ["c","d","e"],
                     [["f","g",["h"],["i",["j","k"],"l"]],"m"],
                     "n",
                     ["o"],
                     [["p"]],
                     [["q"]],
                     [["r","s"]],
                     ["t"]
                    ]
        cls.list1 = list("abcdefghijklmnopqrst")

        # list containing tuples, which will not be flattened
        cls.tree2 = ["a",
                     "b",
                     ("c","d","e"),
                     [["f","g",["h"],["i",["j","k"],"l"]],"m"],
                     "n",
                     ["o"],
                     [["p"]],
                     [["q"]],
                     (["r","s"],),
                     ["t"]
                    ]
        cls.list2 = ["a","b",
                     ("c","d","e"),
                     "f","g","h","i","j","k","l","m","n","o","p","q",
                     (["r","s"],),
                     "t"
                    ] 
 
    def test_parse_list(self):
        # homogeneous lists

        # positive ints
        self.assertEquals(parse_list(str(list(range(50)))),list(range(50)))

        # positive and negative ints
        self.assertEquals(parse_list(str(list(range(-50,100)))),list(range(-50,100)))

        # string
        a = ["asdfads","basfdads","asdfaj3r23fjd8c","adsfj3r92h"]
        self.assertEquals(parse_list(str(a)),a)

        # positive and negative floats
        b = list(numpy.linspace(-5,200,43))
        self.assertEquals(parse_list(str(b)),b)

        # heterogeneous list
        b[10] = True
        b[20] = 3
        b[30] = False
        b[40] = ("this is a test")
        self.assertEquals(parse_list(str(b)),b)

    def test_flatten(self):
        for item1,item2 in zip(flatten_nested_lists_to_generator(self.tree1),self.list1):
            self.assertEquals(item1,item2)

        for item1,item2 in zip(flatten_nested_lists_to_generator(self.list1),self.list1):
            self.assertEquals(item1,item2)

        self.assertNotEquals(self.tree1,self.list1)

        for item1,item2 in zip(flatten_nested_lists_to_generator(self.tree2),self.list2):
            self.assertEquals(item1,item2)

        for item1,item2 in zip(flatten_nested_lists_to_generator(self.list2),self.list2):
            self.assertEquals(item1,item2)

        self.assertNotEquals(self.tree2,self.list2)

    def test_flatten_typesafe(self):
        self.assertRaises(AssertionError,next,flatten_nested_lists_to_generator(5))
        self.assertRaises(AssertionError,next,flatten_nested_lists_to_generator("some_string"))
        self.assertRaises(AssertionError,next,flatten_nested_lists_to_generator(("a","b","c")))

    def test_flatten_nested_lists(self):
        self.assertEquals(flatten_nested_lists_to_list(self.tree1),self.list1)
        self.assertEquals(flatten_nested_lists_to_list(self.list1),self.list1)
        self.assertNotEquals(self.tree1,self.list1)

        self.assertEquals(flatten_nested_lists_to_list(self.tree2),self.list2)
        self.assertEquals(flatten_nested_lists_to_list(self.list2),self.list2)
        self.assertNotEquals(self.tree2,self.list2)

    def test_flatten_nested_lists_typesafe(self):
        self.assertRaises(AssertionError,flatten_nested_lists_to_list,5)
        self.assertRaises(AssertionError,flatten_nested_lists_to_list,"some_string")
        self.assertRaises(AssertionError,flatten_nested_lists_to_list,("a","b","c"))
