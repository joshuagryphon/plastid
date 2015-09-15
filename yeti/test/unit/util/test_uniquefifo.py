#!/usr/bin/env python
"""Tests capabilities of |UniqueFIFO|"""
import unittest
from nose.plugins.attrib import attr

from plastid.util.unique_fifo import UniqueFIFO
from collections import Counter

@attr(test="unit")
class TestUniqueFIFO(unittest.TestCase):

    def test_getitem(self):
        my_fifo = UniqueFIFO(5)
        for i in range(5):
            my_fifo.append(i)
            self.assertEqual(i,my_fifo[-1])
        
        for i in range(5):
            self.assertEqual(i,my_fifo[i])
        
        for i in range(10):
            my_fifo.append(i)
            self.assertEqual(i,my_fifo[-1])
    
    def test_contains(self):
        my_fifo = UniqueFIFO(5)
        self.assertFalse(0 in my_fifo)
        my_fifo.append(0)
        self.assertTrue(0 in my_fifo)
        my_fifo.append(1)
        self.assertTrue(0 in my_fifo)
        my_fifo.append(2)
        self.assertTrue(0 in my_fifo)
        my_fifo.append(3)
        self.assertTrue(0 in my_fifo)
        my_fifo.append(4)
        self.assertTrue(0 in my_fifo)
        my_fifo.append(5)
        self.assertFalse(0 in my_fifo)
    
    def test_iter(self):
        fifo_items1 = list(range(5,10))
        fifo_items2 = list(range(10,15))
        all_fifo_items = fifo_items1 + fifo_items2
        my_fifo = UniqueFIFO(5)
        # add up to length, and make sure at each addition
        # the right element is in each spot
        for item in fifo_items1:
            my_fifo.append(item)
            for c, item in enumerate(my_fifo):
                self.assertEqual(item,fifo_items1[c])
        
        # add at sizes >= max_length, and make sure at each addition
        # the right element is in each spot
        for n,item in enumerate(fifo_items2):
            my_fifo.append(item)
            for c, item in enumerate(my_fifo):
                self.assertEqual(item,all_fifo_items[n+c+1])

    def test_len(self):
        my_fifo = UniqueFIFO(5)
        for i in range(5):
            self.assertEqual(len(my_fifo),i)
            my_fifo.append(i)
        for i in range(5):
            self.assertEqual(len(my_fifo),5)
        
                        
    def test_append_no_exceed_maxlen(self):
        for max_size in range(1,10):
            my_fifo = UniqueFIFO(max_size)
            for i in range(100):
                my_fifo.append(i)
                self.assertTrue(len(my_fifo) <= max_size)
                self.assertTrue(len(list(my_fifo)) <= max_size)

    def test_append_up_to_maxlen(self):
        my_fifo = UniqueFIFO(5)
        for i in range(5):
            my_fifo.append(i)
            self.assertEqual(len(my_fifo),i+1)
    
    def test_append_reorder(self):
        my_fifo = UniqueFIFO(5)
        my_fifo.append(0)
        my_fifo.append(1)
        my_fifo.append(2)
        my_fifo.append(3)
        my_fifo.append(4)
        self.assertListEqual(list(my_fifo),[0,1,2,3,4])
        
        my_fifo.append(3)
        self.assertListEqual(list(my_fifo),[0,1,2,4,3])

        my_fifo.append(0)
        self.assertListEqual(list(my_fifo),[1,2,4,3,0])

        my_fifo.append(4)
        self.assertListEqual(list(my_fifo),[1,2,3,0,4])
        
        my_fifo.append(2)
        self.assertListEqual(list(my_fifo),[1,3,0,4,2])

        my_fifo.append(2)
        self.assertListEqual(list(my_fifo),[1,3,0,4,2])

        my_fifo.append(2)
        self.assertListEqual(list(my_fifo),[1,3,0,4,2])

        my_fifo.append(2)
        self.assertListEqual(list(my_fifo),[1,3,0,4,2])


    def test_append_no_duplicates(self):
        my_fifo = UniqueFIFO(5)
        for i in range(5):
            my_fifo.append(i)
        self.assertListEqual(list(set(my_fifo)),list(my_fifo))
        for i in range(20):
            my_fifo.append(i)
            for v in Counter(my_fifo).values():
                self.assertEqual(v,1)