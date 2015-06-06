#!/usr/bin/env python
"""Test suite for :py:mod:`yeti.util.scriptlib.help_formatters`"""
import unittest
from nose.plugins.attrib import attr
from yeti.util.scriptlib.help_formatters import shorten_help,\
                                                           format_module_docstring,\
                                                           pyrst_pattern,\
                                                           subst_pattern,\
                                                           _separator
@attr(test="unit")
class TestHelpFormatters(unittest.TestCase):
    
    def test_pyrst_pattern_matches(self):
        tests = [(":py:class:`yeti.genomics.genome_array.GenomeArray`",
                  "py",
                  "class",
                  "yeti.genomics.genome_array.GenomeArray"),
                 (":py:mod:`yeti.genomics.genome_array`",
                  "py",
                  "mod",
                  "yeti.genomics.genome_array"),
                (":py:obj:`some_obj`",
                  "py",
                 "obj",
                 "some_obj"),
                 (":class:`yeti.genomics.genome_array.GenomeArray`",
                  None,
                  "class",
                  "yeti.genomics.genome_array.GenomeArray"),
                 (":mod:`yeti.genomics.genome_array`",
                  None,
                  "mod",
                  "yeti.genomics.genome_array"),
                (":obj:`some_obj`",
                 None,
                 "obj",
                 "some_obj"),
                 ]
        for test, domain, role, argument in tests:
            groupdict = pyrst_pattern.search(test).groupdict()
            self.assertEquals(domain,groupdict["domain"])
            self.assertEquals(role,groupdict["role"])
            self.assertEquals(argument,groupdict["argument"])
    
    def test_pyrst_pattern_nonmatches(self):
        tests = [(":py:"),
                 (":py:class:"),
                 ("`genome_array`"),
                 # missing one or more backticks
                 (":py:class:`yeti.genomics.genome_array.GenomeArray"),
                 (":py:class:yeti.genomics.genome_array.GenomeArray`"),
                 (":py:class:yeti.genomics.genome_array.GenomeArray"),
                 # malformed starting token
                 ("py:class:`yeti.genomics.genome_array.GenomeArray`"),
                 ]
        for test in tests:
            self.assertTrue(pyrst_pattern.search(test) is None,
                            "rst pattern matches %s" % test)
    
    def test_subst_pattern_matches(self):
        self.assertEquals(subst_pattern.search("|test_replacement|").groups()[0],"test_replacement")
    
    def test_subst_pattern_nonmatches(self):
        tests = ["|bad_test",
                 "bad_test|"]
        for test in tests:
            self.assertTrue(subst_pattern.search(test) is None,
                            "substitution pattern matches %s" % test)

    def test_shorten_help(self):
        self.assertNotEqual(raw_help1,shortened_help1)
        self.assertEqual(shorten_help(raw_help1),shortened_help1)
    
    def test_format_module_docstring(self):
        expected = _separator + shortened_help1 + _separator
        self.assertEqual(format_module_docstring(raw_help1),expected)
        self.assertEqual(format_module_docstring(shortened_help1),expected)

#===============================================================================
# INDEX: test data
#===============================================================================

raw_help1 = """Implements array-like data structures called |GenomeArray|, that map numerical values 
to genomic positions (nucleotides) under various configurable mapping rules.
Several different implementations are provided, depending on how the alignment
and/or count data is stored:

Important classes
-----------------
    
    |AbstractGenomeArray|
        Base class for all genome arrays, array-like data structures mapping
        read counts to specific nucleotide positions.  
        
        Defines interfaces for:
            - retrieving vectors (as :py:class:`numpy.ndarray` s) of read counts
              at each position of a |GenomicSegment|
              
            - getting the sum of read counts in a dataset
            
            - toggling normalization of fetched counts to reads-per-million
    
    |MutableAbstractGenomeArray|
        Base class for |GenomeArray| and |SparseGenomeArray|. Contains all interfaces
        from |AbstractGenomeArray|, and additionally defines interfaces for:
        
            - setting values, manually or mathematically, over regions of the
              genome, or over the entire genome, element wise
            
See also
--------
:py:meth:`some_method`
    A method that is useful

:py:class:`yeti.genomics.genome_array.GenomeArray`
    A useful class
            
Parameters
----------
foo1 : footype
    Some parameter

foo2 : footype
    another parameter

Raises
------
FooException
    If the foo doesn't bar
"""

shortened_help1 = """Implements array-like data structures called GenomeArray, that map numerical values 
to genomic positions (nucleotides) under various configurable mapping rules.
Several different implementations are provided, depending on how the alignment
and/or count data is stored:

Important classes
-----------------
    
    AbstractGenomeArray
        Base class for all genome arrays, array-like data structures mapping
        read counts to specific nucleotide positions.  
        
        Defines interfaces for:
            - retrieving vectors (as numpy.ndarray s) of read counts
              at each position of a GenomicSegment
              
            - getting the sum of read counts in a dataset
            
            - toggling normalization of fetched counts to reads-per-million
    
    MutableAbstractGenomeArray
        Base class for GenomeArray and SparseGenomeArray. Contains all interfaces
        from AbstractGenomeArray, and additionally defines interfaces for:
        
            - setting values, manually or mathematically, over regions of the
              genome, or over the entire genome, element wise
            
See also
--------
some_method
    A method that is useful

yeti.genomics.genome_array.GenomeArray
    A useful class
"""