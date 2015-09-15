#!/usr/bin/env python
"""Test suite for filters in :py:mod:`plastid.util.io.filters`
"""
import unittest
from nose.plugins.attrib import attr
from plastid.util.services.mini2to3 import cStringIO
from plastid.util.io.filters import SkipBlankReader, CommentReader, FunctionReader


#===============================================================================
# INDEX: test cases
#===============================================================================

@attr(test="unit")
class TestSkipBlankReader(unittest.TestCase):
    """TestCase for :py:class:`SkipBlankReader`"""

    def get_stream(self):
        return SkipBlankReader(cStringIO.StringIO(_NL_TEXT))

    def test_read(self):
        fin = self.get_stream()
        lines = fin.read().strip().split("\n")
        self.assertEqual(len(lines),30)
        for line in lines:
            self.assertNotEqual(line.strip(),"")
        fin.close()

    def test_readlines(self):
        lines = self.get_stream().readlines()
        self.assertEqual(len(lines),30)
        for line in lines:
            self.assertNotEqual(line.strip(),"")

    def test_iter(self):
        for n,line in enumerate(self.get_stream()):
            self.assertNotEqual(line.strip(),"")

        self.assertEqual(n,29)


@attr(test="unit")
class TestCommentReader(TestSkipBlankReader):
    """TestCase for :py:class:`CommentReader`"""
   
    def get_stream(self):
        return CommentReader(cStringIO.StringIO(_COMMENT_TEXT))
    
    def test_get_comments(self):
        reader = self.get_stream()
        comments = reader.get_comments()
        # comments should not be populated until we have read
        self.assertEquals(len(comments),0)
        reader.read()
        # comments should be found afterwards
        self.assertEquals(len(comments),5)
        for comment in comments:
            self.assertTrue(comment.startswith("#"))

@attr(test="unit")
class TestFunctionReader(unittest.TestCase):
    def test_via_backwards(self):
        reader = FunctionReader(cStringIO.StringIO(_NL_TEXT),lambda x: x[::-1])
        for line1, line2 in zip(reader,cStringIO.StringIO(_NL_TEXT)):
            self.assertEqual(line1,line2[::-1])

    def test_via_upper(self):
        reader = FunctionReader(cStringIO.StringIO(_NL_TEXT),str.upper)
        for line1, line2 in zip(reader,cStringIO.StringIO(_NL_TEXT)):
            self.assertEqual(line1,line2.upper())

#===============================================================================
# INDEX: test data
#===============================================================================

_NL_TEXT = """1	Lorem ipsum dolor sit amet, consectetur adipiscing elit. Morbi vitae ipsum vel
     2	nisi dapibus dapibus in vel neque. Proin nec consectetur arcu. Praesent
     3	convallis venenatis metus quis elementum. Ut elementum justo ac faucibus
     4	efficitur. Nam at ornare elit. Sed pulvinar mi sapien, sed faucibus risus
     5	cursus at. Nulla sit amet posuere ex, sit amet convallis ante. Integer
     6	imperdiet purus nec ante pretium fringilla. Vivamus a ligula tristique,
     7	sodales elit et, fringilla ligula. Nulla facilisi. Nunc sed enim non ligula
     8	commodo blandit. Sed ultricies quis urna vel imperdiet. Curabitur leo nisi,
     9	faucibus sed mauris ut, bibendum pretium nibh.
       
    10	Cras libero quam, scelerisque ut suscipit ac, pellentesque non sapien. Vivamus
    11	eu tristique nisi. Sed luctus molestie mollis. Praesent dapibus tincidunt
    12	pretium. Sed faucibus vestibulum est, ac mollis libero dapibus in. Phasellus
    13	nec euismod sapien. Donec faucibus orci sem, vitae hendrerit orci sodales non.
    14	Sed facilisis erat at erat semper facilisis. Proin at orci et ligula mattis
    15	condimentum. Aenean tincidunt, nunc sit amet malesuada scelerisque, nisl augue
    16	imperdiet lectus, et sollicitudin quam quam nec massa. Morbi volutpat nulla et
    17	erat porta tempus.
       
    18	Nulla placerat ipsum elit, at vestibulum urna pellentesque a. Nunc bibendum
    19	convallis orci, vel euismod tellus commodo ullamcorper. Vestibulum vestibulum
    20	lobortis tempor. Maecenas non nisi aliquet lorem tincidunt varius sed aliquam
    21	est. Vivamus egestas, nisi sed laoreet interdum, nisi tellus egestas erat, sit
    22	amet ultricies metus orci sed orci. Aenean at lectus viverra, venenatis mauris
    23	sed, ultricies augue. Sed ac purus vitae dui tempus pulvinar. Duis tincidunt
                               
    24	nisi sit amet purus laoreet placerat. Quisque nulla orci, vestibulum nec metus
    25	eu, convallis consequat libero. Maecenas sed ultrices tellus, non commodo
    26	lectus. Suspendisse non diam eget purus imperdiet aliquam a sed libero.
    									
    27	Quisque eu venenatis tellus. Aliquam ornare purus non ante imperdiet interdum.
    28	Cras pretium, diam ac sagittis luctus, diam metus scelerisque libero, sit amet
                           29	cursus erat nulla ut augue. Ut est tortor, faucibus non nisi a, gravida
    30	fermentum tortor.



"""



_COMMENT_TEXT = """1	Lorem ipsum dolor sit amet, consectetur adipiscing elit. Morbi vitae ipsum vel
     2	nisi dapibus dapibus in vel neque. Proin nec consectetur arcu. Praesent
     3	convallis venenatis metus quis elementum. Ut elementum justo ac faucibus
     4	efficitur. Nam at ornare elit. Sed pulvinar mi sapien, sed faucibus risus
     5	cursus at. Nulla sit amet posuere ex, sit amet convallis ante. Integer
     6	imperdiet purus nec ante pretium fringilla. Vivamus a ligula tristique,
     7	sodales elit et, fringilla ligula. Nulla facilisi. Nunc sed enim non ligula
     8	commodo blandit. Sed ultricies quis urna vel imperdiet. Curabitur leo nisi,
     9	faucibus sed mauris ut, bibendum pretium nibh.
# Here is comment #1  
    10	Cras libero quam, scelerisque ut suscipit ac, pellentesque non sapien. Vivamus
    11	eu tristique nisi. Sed luctus molestie mollis. Praesent dapibus tincidunt
    12	pretium. Sed faucibus vestibulum est, ac mollis libero dapibus in. Phasellus
    13	nec euismod sapien. Donec faucibus orci sem, vitae hendrerit orci sodales non.
    14	Sed facilisis erat at erat semper facilisis. Proin at orci et ligula mattis
    15	condimentum. Aenean tincidunt, nunc sit amet malesuada scelerisque, nisl augue
    16	imperdiet lectus, et sollicitudin quam quam nec massa. Morbi volutpat nulla et
    17	erat porta tempus.
# This is comment #2
    18	Nulla placerat ipsum elit, at vestibulum urna pellentesque a. Nunc bibendum
    19	convallis orci, vel euismod tellus commodo ullamcorper. Vestibulum vestibulum
    20	lobortis tempor. Maecenas non nisi aliquet lorem tincidunt varius sed aliquam
    21	est. Vivamus egestas, nisi sed laoreet interdum, nisi tellus egestas erat, sit
    22	amet ultricies metus orci sed orci. Aenean at lectus viverra, venenatis mauris
    23	sed, ultricies augue. Sed ac purus vitae dui tempus pulvinar. Duis tincidunt
# This is comment #3                               
    24	nisi sit amet purus laoreet placerat. Quisque nulla orci, vestibulum nec metus
    25	eu, convallis consequat libero. Maecenas sed ultrices tellus, non commodo  # this is an end-of-line comment, which should not be removed
    26	lectus. Suspendisse non diam eget purus imperdiet aliquam a sed libero.
# This is comment #4						
    27	Quisque eu venenatis tellus. Aliquam ornare purus non ante imperdiet interdum.
    28	Cras pretium, diam ac sagittis luctus, diam metus scelerisque libero, sit amet
                           29	cursus erat nulla ut augue. Ut est tortor, faucibus non nisi a, gravida
    30	fermentum tortor.
# This is comment #5"""
