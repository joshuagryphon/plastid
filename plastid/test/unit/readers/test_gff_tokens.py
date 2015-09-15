#!/usr/bin/env python
import unittest
import warnings
import copy
from nose.tools import assert_equal, assert_not_equal
from nose.plugins.attrib import attr
from plastid.readers.gff_tokens import escape_GFF3, unescape_GFF3,\
                                                    escape_GTF2, unescape_GTF2,\
                                                    make_GFF3_tokens,\
                                                    make_GTF2_tokens,\
                                                    parse_GFF3_tokens,\
                                                    parse_GTF2_tokens
from plastid.util.services.decorators import catch_warnings

#===============================================================================
# INDEX: test data
#===============================================================================

# string requiring escaping
lipsum = """Lorem %;%ipsum dolor "sit amet, consectetur adipiscing elit, sed do
eiusmod temporincididunt ut labore et dolore magna aliqua. Ut enim ad minim
veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
!*&)$#@(*RVLKJE@*(!U_)$(@*!)JFKJSF)@!(#_!@!(*$9087325gasdfjk@&!*Y*(_::{|}@{!';@
consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum;
dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident,
sunt in culpa qui officia deserunt mollit anim id est laborum.@(#*@FDSFJLK@!="""

lipsum_gff_escaped = """Lorem %25%3B%25ipsum dolor "sit amet%2C consectetur adipiscing elit%2C sed do
eiusmod temporincididunt ut labore et dolore magna aliqua. Ut enim ad minim
veniam%2C quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
!*%26)$#@(*RVLKJE@*(!U_)$(@*!)JFKJSF)@!(#_!@!(*$9087325gasdfjk@%26!*Y*(_::{|}@{!'%3B@
consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum%3B
dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident%2C
sunt in culpa qui officia deserunt mollit anim id est laborum.@(#*@FDSFJLK@!%3D""".replace("\n","%0A")

lipsum_gtf_escaped = """Lorem %25%3B%25ipsum dolor %22sit amet%2C consectetur adipiscing elit%2C sed do
eiusmod temporincididunt ut labore et dolore magna aliqua. Ut enim ad minim
veniam%2C quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
!*%26)$#@(*RVLKJE@*(!U_)$(@*!)JFKJSF)@!(#_!@!(*$9087325gasdfjk@%26!*Y*(_::{|}@{!'%3B@
consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum%3B
dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident%2C
sunt in culpa qui officia deserunt mollit anim id est laborum.@(#*@FDSFJLK@!%3D""".replace("\n","%0A")


# string that escaping should not change
medstring = "A medium sized string with spaces et c"
medstring_escaped = "A medium sized string with spaces et c"

# tokens formatted in GFF3

_GFF3_REPEATED = [
                ("repeated token no trailing semicolon"    ,'a=1;c=3;b=2;e=5;d=4;z=26;c=7'),
                ("repeated token trailing semicolon"       ,'a=1;c=3;b=2;e=5;d=4;z=26;c=7;'),
                  ]

_GTF2_REPEATED = [
                ("repeated token no trailing semicolon"    ,'a "1"; c "3"; b "2"; e "5"; d "4"; z "26"; c "7"'),
                ("repeated token trailing semicolon"       ,'a "1"; c "3"; b "2"; e "5"; d "4"; z "26"; c "7";'),
                  ]

_GFF3_TOKENS = [("simple values no trailing semicolon"     ,'a=1;c=3;b=2;e=5;d=4;z=26'),
                ("too much spacing no trailing semicolon"  ,'a=1;   c=3;b=2  ;e=5;d=4;z=26'),
                ("escaped strings no trailing semicolon"   ,'a=5;note1=%s;note2=%s;note3=%s' % (medstring,medstring_escaped,lipsum_gff_escaped)),
                ("simple values with trailing semicolon"   ,'a=1;c=3;b=2;e=5;d=4;z=26;'),
                ("too much spacing with trailing semicolon",'a=1;   c=3;b=2  ;e=5;d=4;z=26;'),
                ("two value non-expanded"                  ,'a=1;c=3,7;b=2;e=5;d=4;z=26;'),
                ("escaped strings with trailing semicolon" ,'a=5;note1=%s;note2=%s;note3=%s;' % (medstring,medstring_escaped,lipsum_gff_escaped)),
                ("escaped key"                             ,'c%%2Cd=5;a=5;note1=%s;note2=%s;note3=%s;' % (medstring,medstring_escaped,lipsum_gff_escaped)),
                
                # keep GFF3-specific elements at the end, so they are not tripped in GTF2
                ("two value GFF3 expanded no trailing semicolon",'Parent=p1,p2;a=1;c=3,7;b=2;e=5;d=4;z=26;Alias=alias1,alias2'),
                ("two value GFF3 expanded with trailing semicolon",'Parent=p1,p2;a=1;c=3,7;b=2;e=5;d=4;z=26;Alias=alias1,alias2;'),
                ]

# tokens formatted in GTF2
_GTF2_TOKENS = [("simple values no trailing semicolon"     ,'a "1"; c "3"; b "2"; e "5"; d "4"; z "26"'),
                ("too much spacing no trailing semicolon"  ,'a     "1"; c    "3"; b "2"; e   "5"; d "4"; z "26"'),
                ("escaped strings no trailing semicolon"   ,'a "5"; note1 "%s"; note2 "%s"; note3 "%s"' % (medstring,medstring_escaped,lipsum_gtf_escaped)),
                ("simple values with trailing semicolon"   ,'a "1"; c "3"; b "2"; e "5"; d "4"; z "26;"'),
                ("too much spacing with trailing semicolon",'a     "1"; c    "3"; b "2"; e   "5"; d "4"; z "26";'),
                ("two value non-expanded"                  ,'a "1"; c "3,7"; b "2"; e "5"; d "4"; z "26";'),
                ("escaped strings with trailing semicolon" ,'a "5"; note1 "%s"; note2 "%s"; note3 "%s";' % (medstring,medstring_escaped,lipsum_gtf_escaped)),
                ("escaped key"                             ,'c%%2Cd "5"; a "5"; note1 "%s"; note2 "%s"; note3 "%s";' % (medstring,medstring_escaped,lipsum_gtf_escaped)), # key must be escaped
                ]

_TOKEN_PARSE_REPEATED = [
                dict(a="1",c=["3","7"],b="2",e="5",d='4',z="26"), # repeated tokens, must be catenated
                dict(a="1",c=["3","7"],b="2",e="5",d='4',z="26"), # repeated tokens, must be catenated
                         ]

_TOKEN_WRITE_REPEATED = [
                dict(a="1",c="3,7",b="2",e="5",d='4',z="26"), # repeated tokens, must be catenated
                dict(a="1",c="3,7",b="2",e="5",d='4',z="26"), # repeated tokens, must be catenated
                         ]

# dicitonaries corresponding to token strings above
_GFF3_TOKEN_DICTS = [dict(a="1",c="3",b="2",e="5",d='4',z="26"), # just values, no trailing semicolon
                dict(a="1",c="3",b="2",e="5",d='4',z="26"), # too much spacing
                dict(a="5",note1=medstring,                 # requires escaping
                     note2=medstring,
                     note3=lipsum),
                dict(a="1",c="3",b="2",e="5",d='4',z="26"), # just values, trainling semicolon
                dict(a="1",c="3",b="2",e="5",d='4',z="26"), # too much spacing
                dict(a="1",c="3,7",b="2",e="5",d='4',z="26"), # two value non-expanded
                dict(a="5",note1=medstring,                 # requires escaping
                     note2=medstring,
                     note3=lipsum),
                {  "c,d"   : "5",
                   "a"     : "5",
                   "note1" : medstring,
                   "note2" : medstring,
                   "note3" : lipsum
                 },
                dict(Parent=["p1","p2"],
                     Alias=["alias1","alias2"],
                     a="1",c="3,7",b="2",e="5",d='4',z="26"),
                dict(Parent=["p1","p2"],
                     Alias=["alias1","alias2"],
                     a="1",c="3,7",b="2",e="5",d='4',z="26"),
                ]

_GTF2_TOKEN_DICTS = _GFF3_TOKEN_DICTS[:-2]

#===============================================================================
# INDEX: tests for token parsing & writing
#===============================================================================

@attr(test="unit")
def test_escape_GFF3():
    yield assert_equal, escape_GFF3(lipsum), lipsum_gff_escaped
    yield assert_equal, escape_GFF3(medstring), medstring_escaped
    yield assert_not_equal, lipsum, lipsum_gff_escaped

@attr(test="unit")
def test_unescape_GFF3():
    yield assert_equal, unescape_GFF3(lipsum_gff_escaped), lipsum
    yield assert_equal, unescape_GFF3(medstring_escaped), medstring
    
@attr(test="unit")
def test_escape_GTF2():
    yield assert_equal, escape_GTF2(lipsum), lipsum_gtf_escaped
    yield assert_equal, escape_GTF2(medstring), medstring_escaped
    yield assert_not_equal, lipsum, lipsum_gtf_escaped

@attr(test="unit")
def test_unescape_GTF2():
    yield assert_equal, unescape_GTF2(lipsum_gtf_escaped), lipsum
    yield assert_equal, unescape_GTF2(medstring_escaped), medstring


@attr(test="unit")
class TestGFF3_TokenParsing(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.fmt    = "GFF3"
        cls.tokens = _GFF3_TOKENS
        cls.token_dicts = _GFF3_TOKEN_DICTS
        cls.repeated = _GFF3_REPEATED
        cls.tokenizer = staticmethod(make_GFF3_tokens)
        cls.parser = staticmethod(parse_GFF3_tokens)
        
    def test_token_parse(self):
        pass
        """Test parsing of tokens in *attributes* column of GFF"""
        # test that all elements are in dict
        for c, (my_dict, (my_title,my_tokens)) in enumerate(zip(self.token_dicts,self.tokens)):
            found = self.parser(my_tokens)
            msg = "Failed %s parsing test '%s'\n    Expected: %s\n    Found: %s" % (self.fmt,my_title,my_dict,found)
            self.assertEqual(my_dict,found,msg)
         
        self.assertEqual(c,len(self.tokens)-1,"Not all tokens parsed during test. Expected %s. Parsed %s" % (len(self.tokens)-1,c))
 
    def test_token_write_all(self,dicts=None):
        pass
        """Test export of GFF3 tokens
         
        Parameters
        ----------
        dicts : list<dict>
            List of dictionaries to convert to GFF3 tokens
        """
        if dicts is None:
            dicts = self.token_dicts
            
        # test total export from written tokens
        for c, my_dict in enumerate(dicts):
            my_str = self.tokenizer(my_dict,escape=True)
            new_dict = self.parser(my_str)
            self.assertEqual(my_dict,new_dict)
         
        self.assertGreater(c,0)
        
    def test_token_write_subset(self,dicts=None,non_excludes=[]):
        pass
        """Test export of GFF3 tokens, excluding random subsets
        
        Parameters
        ----------
        dicts : list<dict>
            List of dictionaries to convert to GFF3 tokens
        
        non_excludes : list<str>
            Keys that cannot be excluded from export
        """        
        # test excluding random key subsets from written tokens
        if dicts is None:
            dicts = self.token_dicts
            
        for c, old_dict in enumerate(dicts):

            excludes = set(list(old_dict.keys())[0:3]) - set(non_excludes)
            my_str = self.tokenizer(old_dict,excludes=list(excludes),escape=True)
            new_dict = self.parser(my_str)
            
            old_keys = set(old_dict.keys())
            new_keys = set(new_dict.keys())
            self.assertEqual(old_keys-excludes,new_keys)
            for k in (new_keys):
                self.assertEqual(old_dict[k],new_dict[k])
            
        self.assertGreater(c,0)

    def test_repeated_token_raises_warning(self):
        with warnings.catch_warnings(record=True) as warns:
             for c,(name, inp) in enumerate(self.repeated):
                warnings.simplefilter("always")
                msg = "%s parser failed to raise warning for %s" % (self.fmt, name)
                self.parser(inp)
                self.assertEqual(len(warns),c+1)

        
@attr(test="unit")
class TestGTF2_TokenParsing(TestGFF3_TokenParsing):
    """Test GTF2 feature reading and writing.
  
    See http://mblab.wustl.edu/GTF22.html
    """
  
    @classmethod
    def setUpClass(cls):
        cls.fmt    = "GTF2"
        cls.tokens = _GTF2_TOKENS
        cls.token_dicts = _GTF2_TOKEN_DICTS
        cls.repeated = _GTF2_REPEATED
        cls.tokenizer = staticmethod(make_GTF2_tokens)
        cls.parser = staticmethod(parse_GTF2_tokens)
     
    @staticmethod
    def modify_dict(my_dict):
        """Add *gene_id* and *transcript_id* attributes to attribute dictionary,
        because these are required by GTF2 elements
         
        Parameters
        ----------
        my_dict : dict
            Key-value pairs to export
         
        Returns
        -------
        dict
        """
        my_dict = copy.deepcopy(my_dict)
        my_dict["gene_id"] = "some_gene"
        my_dict["transcript_id"] = "some_transcript"
        return my_dict
    
    def test_token_write_all(self):
        my_dicts = [self.modify_dict(X) for X in self.token_dicts]
        super(TestGTF2_TokenParsing,self).test_token_write_all(dicts=my_dicts)
             
    def test_token_write_subset(self):
        my_dicts = [self.modify_dict(X) for X in self.token_dicts]
        super(TestGTF2_TokenParsing,self).test_token_write_subset(dicts=my_dicts,non_excludes=["gene_id","transcript_id"])

    @catch_warnings("ignore") # need to catch warnings because duplicate GTF2 tags raise FileFormatWarnings
    def test_gtf2_list_duplicate_tags(self):
        test_cases = [('gene_id "mygene"; transcript_id "mytranscript"; tag "tag value";',
                       {'gene_id' : 'mygene', 'tag' : 'tag value', 'transcript_id' : 'mytranscript'}),
                      ('gene_id "mygene"; transcript_id "mytranscript"; tag "tag value"; tag "tag value 2";',
                       {'gene_id' : 'mygene', 'tag' : 'tag value,tag value 2', 'transcript_id' : 'mytranscript'}),
                       ('gene_id "mygene"; tag "tag value"; transcript_id "mytranscript"; tag "tag value 2";',
                       {'gene_id' : 'mygene', 'tag' : 'tag value,tag value 2', 'transcript_id' : 'mytranscript'}),
                     ]
        for inp, outp in test_cases:
            self.assertDictEqual(self.parser(inp),outp)
