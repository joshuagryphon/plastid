#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.util.services.sets`"""

from plastid.util.services.sets import merge_sets
from nose.plugins.attrib import attr
from nose.tools import assert_set_equal

def check_merge_equal(expected, found):
    yield assert_set_equal, expected, found

@attr(test="unit")
def test_merge_sets(fn=merge_sets):
    """Tests :py:func:`merge_sets`"""
    a = [set(['h']),
         set(['a', 'b', 'm']),
         set(['c']),
         set(['c']),
         set(['b', 'j']),
         set(['i']),
         set(['k', 'o']),
         set(['n']),
         set(['g', 'o']),
         set(['i', 'k']),
         set(['e', 'i']),
         set(['a']),
         set(['d', 'h']),
         set(['l']),
         set(['g', 'j', 'k']),
         set(['f']),
         set(['b']),
         set(['g', 'n']),
         set(['d', 'm', 'p']),
         set(['i', 'n'])]
    b = [set(
            ['g', 'm']),
         set(['o']),
         set(['a', 'g']),
         set(['f']),
         set(['h']),
         set(['h']),
         set(['o']),
         set(['b', 'l']),
         set(['e']),
         set(['p']),
         set(['j']),
         set(['p']),
         set(['k']),
         set(['c', 'f']),
         set(['b', 'c']),
         set(['b']),
         set(['i', 'k']),
         set(['g']),
         set(['j', 'm']),
         set(['b', 'n']),
         set(['i']),
         set(['d', 'o']),
         set(['b', 'l']),
         set(['a']),
         set(['a', 'p'])]
    _tform = lambda x: frozenset([frozenset(X) for X in x])
    yield check_merge_equal, _tform(fn(a)),  _tform([set(['c']),
                     set(['f']),
                     set(['l']),
                     set(['a', 'b', 'd', 'e', 'g', 'h', 'i', 'j', 'k', 'm', 'n', 'o', 'p'])])
    yield check_merge_equal, _tform(fn(b)), _tform([set(['h']),
                     set(['e']),
                     set(['i', 'k']),
                     set(['d', 'o']),
                     set(['b', 'c', 'f', 'l', 'n']),
                     set(['a', 'g', 'j', 'm', 'p'])])  
