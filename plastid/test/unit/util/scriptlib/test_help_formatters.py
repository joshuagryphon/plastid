#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.util.scriptlib.help_formatters`"""
from nose.tools import assert_equal, assert_not_equal, assert_true
from nose.plugins.attrib import attr
from plastid.util.scriptlib.help_formatters import shorten_help,\
                                                           format_module_docstring,\
                                                           pyrst_pattern,\
                                                           subst_pattern,\
                                                           link_pattern,\
                                                           _separator
@attr(test="unit")
class TestHelpFormatters():
    
    def test_pyrst_pattern_matches(self):
        tests = [(":py:class:`plastid.genomics.genome_array.GenomeArray`",
                  "py",
                  "class",
                  "plastid.genomics.genome_array.GenomeArray",
                  None),
                 (":py:mod:`plastid.genomics.genome_array`",
                  "py",
                  "mod",
                  "plastid.genomics.genome_array",
                  None),
                (":py:obj:`some_obj`",
                  "py",
                 "obj",
                 "some_obj",
                  None),
                 (":class:`plastid.genomics.genome_array.GenomeArray`",
                  None,
                  "class",
                  "plastid.genomics.genome_array.GenomeArray",
                  None),
                 (":mod:`plastid.genomics.genome_array`",
                  None,
                  "mod",
                  "plastid.genomics.genome_array",
                  None),
                (":obj:`some_obj`",
                 None,
                 "obj",
                 "some_obj",
                  None),
                (":obj:`some_obj <with pointer>`",
                 None,
                 "obj",
                 "some_obj",
                 "with pointer"),
                (":domain:role:`some_arg <with pointer>`",
                 "domain",
                 "role",
                 "some_arg",
                 "with pointer"),
                 ]
        for test, domain, role, argument, pointer in tests:
            groupdict = pyrst_pattern.search(test).groupdict()
            yield self.check_pyrst_pattern,domain,groupdict["domain"]
            yield self.check_pyrst_pattern,role,groupdict["role"]
            yield self.check_pyrst_pattern,argument,groupdict["argument"]
            yield self.check_pyrst_pattern,pointer,groupdict["pointer"]
    
    @staticmethod
    def check_pyrst_pattern(expected,found,message=""):
        if len(message) == 0:
            message = "Expected %s, found %s." % (expected,found)
        assert_equal(expected,found,message)
        
    def test_pyrst_pattern_nonmatches(self):
        tests = [(":py:"),
                 (":py:class:"),
                 ("`genome_array`"),
                 # missing one or more backticks
                 (":py:class:`plastid.genomics.genome_array.GenomeArray"),
                 (":py:class:plastid.genomics.genome_array.GenomeArray`"),
                 (":py:class:plastid.genomics.genome_array.GenomeArray"),
                 # malformed starting token
                 ("py:class:`plastid.genomics.genome_array.GenomeArray`"),
                 ]
        for test in tests:
            yield self.check_pyrst_pattern, pyrst_pattern.search(test), None, "rst pattern matches %s" % test
                
    def test_subst_pattern_matches(self):
        assert_equal(subst_pattern.search("|test_replacement|").groups()[0],"test_replacement")
    
    def test_subst_pattern_nonmatches(self):
        tests = ["|bad_test",
                 "bad_test|"]
        for test in tests:
            assert_true(subst_pattern.search(test) is None,
                            "substitution pattern matches %s" % test)

    def test_link_pattern_matches(self):
        assert_equal(link_pattern.search("`My link`_").groups()[0],"My link")
        assert_equal(link_pattern.search("`Mylink`_").groups()[0],"Mylink")
        assert_equal(link_pattern.search("`Mylink (with parentheses)`_").groups()[0],"Mylink (with parentheses)")
        assert_equal(link_pattern.search("`Mylink <with_subs>`_").groups()[0],"Mylink")
    
    def test_link_pattern_nonmatches(self):
        tests = ["|bad_test",
                 "bad_test|"]
        for test in tests:
            assert_true(link_pattern.search(test) is None,
                            "substitution pattern matches %s" % test)

    def test_shorten_help(self):
        assert_not_equal, raw_help1,shortened_help1
        assert_equal(shorten_help(raw_help1), shortened_help1)
    
    def test_format_module_docstring(self):
        expected = _separator + "\n" + shortened_help1 + "\n" + _separator
        assert_equal(format_module_docstring(raw_help1), expected)
        assert_equal(format_module_docstring(shortened_help1),expected)

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

:py:class:`plastid.genomics.genome_array.GenomeArray`
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

plastid.genomics.genome_array.GenomeArray
    A useful class
"""