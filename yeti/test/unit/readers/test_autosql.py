#!/usr/bin/env python
"""Test suites for :py:mod:`plastid.readers.autosql`"""
import unittest
import itertools
import nose.tools
from nose.plugins.attrib import attr
from collections import OrderedDict
from plastid.util.services.decorators import skip_if_abstract
from plastid.readers.autosql import AutoSqlField, SizedAutoSqlField, ValuesAutoSqlField, AutoSqlDeclaration, AbstractAutoSqlElement


_autosql_declarations = [
"""table easy_table
"A table with a comment on the next line" 
    (
    uint number auto; "a number with a token"
    uint [3] points ; "r,g,b values"
    lstring  my_string ; "a long string"
    uint a_field_size ; "the size for the next field"
    float [a_field_size] float_array ; "an array of floats"
    set(a,b,c) alpha ; "the first three letters of the alphabet"
    )
""",
"""table easy_table "A table with a comment on the next line" 
    (
    uint number auto; "a number with a token"
    uint [3] points ; "r,g,b values"
    lstring  my_string ; "a long string"
    uint a_field_size ; "the size for the next field"
    float [a_field_size] float_array ; "an array of floats"
    set(a,b,c) alpha ; "the first three letters of the alphabet"
    )
""",
"""table easy_table "A table with a comment on the next line" ( 
uint number auto; "a number with a token"
uint [3] points ; "r,g,b values"
lstring  my_string ; "a long string"
uint a_field_size ; "the size for the next field"
float [a_field_size] float_array ; "an array of floats"
set(a,b,c) alpha ; "the first three letters of the alphabet")
"""
]


match_texts = { "AutoSqlField" : ['uint red ; "red rgb value"',
                                  'uint red;  "red rgb value"',
                                  '    uint red ;  "red rgb value"',
                                  'uint red;  "red rgb value"    ',
                                  'float orange ; "orange value"',
                                  'string  desc    ; "description"',
                                  'lstring  desc    ; "description"',
                                  'float orange primary ; "orange value"',
                                  'float orange primary index ; "orange value"',
                                  'float orange primary index[5] ; "orange value"',
                                  'float orange primary index [5] ; "orange value"',
                                  'float orange primary index [5] auto ; "orange value"',
                                  'float orange primary index[5] auto ; "orange value"',
                                  'float orange primary index [ 5 ] auto ; "orange value"',
                                  'float orange auto primary index [5] ; "orange value"',
                                  'float orange auto primary index; "orange value"',
                                  'float orange primary auto index; "orange value"',
                                  'float orange index primary auto; "orange value"',
                                  'float orange index [ 5] primary auto; "orange value"',
                                 ],
               
                "SizedAutoSqlField"  : ['uint [ 3 ] points ; "r,g,b values"',
                                        'uint [3] points ; "r,g,b values"',
                                        'uint[3] points ; "r,g,b values"',
                                        'uint[3] points; "r,g,b values"',
                                        'float [ 2 ] cmplx ; "real and imaginary parts"',
                                        'string [10] text ; "some text"',
                                        'uint [ numpoints ] points ; "r,g,b values"',
                                        'uint [ 3 ] points primary ; "r,g,b values"',
                                        'uint [ 3 ] points index ; "r,g,b values"',
                                        'uint [ 3 ] points auto ; "r,g,b values"',
                                        'uint [ 3 ] points primary index auto ; "r,g,b values"',
                                        'uint [ 3 ] points index auto primary ; "r,g,b values"',
                                        'uint [ 3 ] points primary index auto; "r,g,b values"',
                                        'uint [ 3 ] points primary auto ; "r,g,b values"',
                                        'uint [ 3 ] points index primary ; "r,g,b values"',
                                        'uint [ 3 ] points primary index [12] auto ; "r,g,b values"',
                                        'uint [ 3 ] points auto primary index [12]; "r,g,b values"',
                                        ],
                
                "ValuesAutoSqlField" : ['set (a,b,c) alphabet ; "the first three letters"',
                                        'set ( a , b , c ) alphabet ; "the first three letters"',
                                        'set ( a , b , c ) alphabet; "the first three letters"',
                                        'set (a, b, c) alphabet ; "the first three letters"',
                                        'enum (a,b,c) alphabet ; "the first three letters"',
                                        'enum ( a , b , c ) alphabet ; "the first three letters"',
                                        'enum ( a , b , c ) alphabet; "the first three letters"',
                                        'enum (a, b, c) alphabet ; "the first three letters"',
                                        'set (a,b,c) alphabet primary ; "the first three letters"',
                                        'set (a,b,c) alphabet index ; "the first three letters"',
                                        'set (a,b,c) alphabet primary index [2] ; "the first three letters"',
                                        'set (a,b,c) alphabet index [2] primary ; "the first three letters"',
                                        ],
                
                "AutoSqlDeclaration" : _autosql_declarations,
               }
"""Text snippets that should be matched by the RegExes in various classes"""

@attr(test="unit")
def test_mask_comments():
    # test comment masking on a variety of fields and on a declaration
    # make sure strings are returned as expected, and make sure location
    # of comment in fields matches index of first quotation mark
    field_tests = ['uint red ; "red rgb value"',
                   'uint red;  "red rgb value"',
                   '    uint red ;  "red rgb value"',
                   'float orange primary index[5] ; "orange value"',
                   'float orange primary index [5] ; "orange value"',
                   'float orange primary index [5] auto ; "orange value"',
                   'float orange primary index[5] auto ; "orange value"',                   
                   'enum ( a , b , c ) alphabet ; "the first three letters"',
                   'enum ( a , b , c ) alphabet; "the first three letters"',
                   'enum (a, b, c) alphabet ; "the first three letters"',
                   'set (a,b,c) alphabet primary ; "the first three letters"',
                   'set (a,b,c) alphabet index ; "the first three letters"',
                   'set (a,b,c) alphabet primary index [2] ; "the first three letters"',                   
                   'uint[3] points ; "r,g,b values"',
                   'uint[3] points; "r,g,b values"',
                   'float [ 2 ] cmplx ; "real and imaginary parts"',
                   'string [10] text ; "some text"',
                   'uint [ numpoints ] points ; "r,g,b values"',
                   'uint [ 3 ] points primary ; "r,g,b values"',
                   'uint [ 3 ] points index ; "r,g,b values"',
                   'uint [ 3 ] points auto ; "r,g,b values"',                  
                   ]
    field_test_results = ['uint red ; "xxxxxxxxxxxxx"',
                   'uint red;  "xxxxxxxxxxxxx"',
                   '    uint red ;  "xxxxxxxxxxxxx"',
                   'float orange primary index[5] ; "xxxxxxxxxxxx"',
                   'float orange primary index [5] ; "xxxxxxxxxxxx"',
                   'float orange primary index [5] auto ; "xxxxxxxxxxxx"',
                   'float orange primary index[5] auto ; "xxxxxxxxxxxx"',                   
                   'enum ( a , b , c ) alphabet ; "xxxxxxxxxxxxxxxxxxxxxxx"',
                   'enum ( a , b , c ) alphabet; "xxxxxxxxxxxxxxxxxxxxxxx"',
                   'enum (a, b, c) alphabet ; "xxxxxxxxxxxxxxxxxxxxxxx"',
                   'set (a,b,c) alphabet primary ; "xxxxxxxxxxxxxxxxxxxxxxx"',
                   'set (a,b,c) alphabet index ; "xxxxxxxxxxxxxxxxxxxxxxx"',
                   'set (a,b,c) alphabet primary index [2] ; "xxxxxxxxxxxxxxxxxxxxxxx"',                   
                   'uint[3] points ; "xxxxxxxxxxxx"',
                   'uint[3] points; "xxxxxxxxxxxx"',
                   'float [ 2 ] cmplx ; "xxxxxxxxxxxxxxxxxxxxxxxx"',
                   'string [10] text ; "xxxxxxxxx"',
                   'uint [ numpoints ] points ; "xxxxxxxxxxxx"',
                   'uint [ 3 ] points primary ; "xxxxxxxxxxxx"',
                   'uint [ 3 ] points index ; "xxxxxxxxxxxx"',
                   'uint [ 3 ] points auto ; "xxxxxxxxxxxx"',                  
                   ]
    
    declaration_test     = 'table my_table "a test table for me " (\n' + "\n".join(field_tests) + "\n)"
    declaration_expected = 'table my_table "xxxxxxxxxxxxxxxxxxxx" (\n' + "\n".join(field_test_results) + "\n)"
    for test, expected in zip(field_tests,field_test_results):
        found = AbstractAutoSqlElement.mask_comments(test)
        expected_start = test.index('"')
        nose.tools.assert_equal(expected,found[0],
                                "Failed to mask field comment '%s'. Expected '%s', got '%s'" % (test,expected,found[0]))
        nose.tools.assert_equal(expected_start+1,found[1][0][0])

    declaration_found = AbstractAutoSqlElement.mask_comments(declaration_test)
    nose.tools.assert_equal(declaration_expected,declaration_found[0],
                            "Failed to mask comment '%s'.Expected:\n'%s'\n Got:\n'%s'" % (declaration_test,
                                                                                      declaration_expected,
                                                                                      declaration_found[0]))
@attr(test="unit")
class TestAbstractAutoSqlElement(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # define all of these in subclasses
        cls.test_class = None
        
        # instances of autoSql text that match cls.test_class.match_pattern
        cls.matching_text     = []
        
        # instances of autoSql text that do not match cls.test_class.match_pattern
        cls.non_matching_text = []
        
        # dictionary mapping indices of autoSql specs in cls.matching_text
        # to a list of tuples of unparsed text and expected result objects
        cls.expected_values = {}

        # record holding values        
        cls.test_record = OrderedDict()

    @skip_if_abstract
    def test_match_pattern(self):
        for c,text in enumerate(self.matching_text):
            self.assertTrue(self.test_class.match_pattern.search(text) is not None,
                            "%s failed to match text '%s'" % (self.test_class.__name__,text))
        
        self.assertGreater(c+1,0,"%s no matches tested in test_match_pattern" % self.test_class.__name__)

    @skip_if_abstract
    def test_not_match_pattern(self):
        for c,text in enumerate(self.non_matching_text):
            self.assertFalse(self.test_class.match_pattern.search(text) is not None,
                             "%s matched text that should have been unmatchable: %s" % (self.test_class.__name__,text))

        self.assertGreater(c+1,0,"%s no matches tested in test_not_match_pattern" % self.test_class.__name__)

    @skip_if_abstract
    def test_add_type(self):
        my_instance = self.test_class(self.matching_text[0])
        fname = "new_type"
        self.assertNotIn(fname,my_instance.field_types)
        my_instance.add_type(fname,int)
        self.assertIn(fname,my_instance.field_types)
        self.assertEqual(my_instance.field_types[fname],int)

    @skip_if_abstract
    def test_call(self):
        for n, declaration in enumerate(self.matching_text):
            parser = self.test_class(declaration)
            for test_str, expected in self.expected_values[n]:
                found = parser(test_str,self.test_record)
                self.assertEqual(found,expected,
                                 "%s parser results incorrect. Got '%s', expected '%s'" % (self.test_class.__name__,
                                                                                            found,
                                                                                            expected))
@attr(test="unit")
class TestAutoSqlField(TestAbstractAutoSqlElement):
    
    @classmethod
    def setUpClass(cls):
        # define all of these in subclasses
        cls.test_class = AutoSqlField
        cls.match_key = cls.test_class.__name__
        
        # instances of autoSql text that match cls.test_class.match_pattern
        cls.matching_text = match_texts[cls.match_key]
        
        # instances of autoSql text that do not match cls.test_class.match_pattern
        cls.non_matching_text = list(itertools.chain(*map(lambda x: x[1],filter(lambda x: x[0]!= cls.match_key,match_texts.items()))))

        cls.test_record = OrderedDict()
        
        # list of lists of tuples of unparsed text and expected result objects
        # for each autoSql declaration in cls.matching_text
        cls.expected_values = [[("3",3),
                                (" 0",0),
                                ("20",20),
                                ],
                               [("3",3),
                                (" 0",0),
                                ("20",20),
                                ],
                               [("3",3),
                                (" 0",0),
                                ("20",20),
                                ],
                               [("3",3),
                                (" 0",0),
                                ("20",20),
                                ],
                               [("3.0",3.0),
                                (" 0.0",0.0),
                                ("-275.4",-275.4)],
                               [("some text with spaces","some text with spaces"),
                                ("some text,with,commas","some text,with,commas"),
                                ("asf32Faf!#%#@%","asf32Faf!#%#@%")],
                               [("some text with spaces","some text with spaces"),
                                ("some text,with,commas","some text,with,commas"),
                                ("asf32Faf!#%#@%","asf32Faf!#%#@%")],
                              ] +  ([[("3.0",3.0),
                                      (" 0.0",0.0),
                                      ("-275.4",-275.4)]]*12)

@attr(test="unit")
class TestSizedAutoSqlField(TestAbstractAutoSqlElement):
    
    @classmethod
    def setUpClass(cls):
        # define all of these in subclasses
        cls.test_class = SizedAutoSqlField
        cls.match_key = cls.test_class.__name__
        
        # instances of autoSql text that match cls.test_class.match_pattern
        cls.matching_text = match_texts[cls.match_key]
        
        # instances of autoSql text that do not match cls.test_class.match_pattern
        cls.non_matching_text = list(itertools.chain(*map(lambda x: x[1],filter(lambda x: x[0]!= cls.match_key,match_texts.items()))))
        
        cls.test_record = OrderedDict(numpoints=3)

        # list of lists of tuples of unparsed text and expected result objects
        # for each autoSql declaration in cls.matching_text
        cls.expected_values = [[("3,4,5",(3,4,5)),
                                ("3, 4, 5",(3,4,5)),
                                ("   201,0, 476",(201,0, 476)),
                                ],
                               [("3,4,5",(3,4,5)),
                                ("3, 4, 5",(3,4,5)),
                                ("   201,0, 476",(201,0, 476)),
                                ],
                               [("3,4,5",(3,4,5)),
                                ("3, 4, 5",(3,4,5)),
                                ("   201,0, 476",(201,0, 476)),
                                ],
                               [("3,4,5",(3,4,5)),
                                ("3, 4, 5",(3,4,5)),
                                ("   201,0, 476",(201,0, 476)),
                                ],
                               [("5.4,3.2",(5.4,3.2)),
                                ("5.4,3.2   ",(5.4,3.2)),
                                ("    5.4,3.2",(5.4,3.2)),
                                ("5.4,   3.2",(5.4,3.2)),
                                ],
                               [("abcdefghij","abcdefghij")],
                               [("3,4,5",(3,4,5)),
                                ("3, 4, 5",(3,4,5)),
                                ("   201,0, 476",(201,0, 476)),
                                ],
                              ] + ([[("3,4,5",(3,4,5)),
                                ("3, 4, 5",(3,4,5)),
                                ("   201,0, 476",(201,0, 476)),
                                ]]*10)
    
@attr(test="unit")
class TestValuesAutoSqlField(TestAbstractAutoSqlElement):

    @classmethod
    def setUpClass(cls):
        # define all of these in subclasses
        cls.test_class = ValuesAutoSqlField
        cls.match_key = cls.test_class.__name__
        
        # instances of autoSql text that match cls.test_class.match_pattern
        cls.matching_text = match_texts[cls.match_key]
        
        # instances of autoSql text that do not match cls.test_class.match_pattern
        cls.non_matching_text = list(itertools.chain(*map(lambda x: x[1],filter(lambda x: x[0]!= cls.match_key,match_texts.items()))))
        
        cls.test_record = OrderedDict(r=255,g=5,c=10)
        
        # list of lists of tuples of unparsed text and expected result objects
        # for each autoSql declaration in cls.matching_text
        cls.expected_values = [[(("a,b,c"),set(list("abc"))),
                                (("a,c,b"),set(list("abc"))),
                                (("a,b"),set(list("ab"))),
                                (("a"),set(list("a"))),
                                ((""),set()),
                                (("    "),set()),
                                (("a, b, c"),set(list("abc"))),
                                ((" a , b , c "),set(list("abc"))),
                                ((" a,b, c"),set(list("abc"))),
                                ],
                              ]*12

@attr(test="unit")
class TestAutoSqlDeclaration(TestAbstractAutoSqlElement):

    @classmethod
    def setUpClass(cls):
        # define all of these in subclasses
        cls.test_class = AutoSqlDeclaration
        cls.match_key = cls.test_class.__name__
        
        # instances of autoSql text that match cls.test_class.match_pattern
        cls.matching_text = match_texts[cls.match_key]
        
        # instances of autoSql text that do not match cls.test_class.match_pattern
        cls.non_matching_text = list(itertools.chain(*map(lambda x: x[1],filter(lambda x: x[0]!= cls.match_key,match_texts.items()))))
        
        cls.test_record = OrderedDict()
        
        # list of lists of tuples of unparsed text and expected result objects
        # for each autoSql declaration in cls.matching_text
        cls.expected_values = [[(("3    1,2,3    my string with spaces    5    1.1,1.2,1.3,1.4,1.5    a,b,c".replace("    ","\t")),
                                 OrderedDict([("number",3),
                                              ("points",(1,2,3)),
                                              ("my_string","my string with spaces"),
                                              ("a_field_size",5),
                                              ("float_array",(1.1,1.2,1.3,1.4,1.5)),
                                              ("alpha",set(list("abc")))])),
                               
                                (("200    5, 212341,12341323    my string with spaces and @&$@(#!@&!    2    3.0,-4586.2    a".replace("    ","\t")),
                                 OrderedDict([("number",200),
                                              ("points",(5,212341,12341323)),
                                              ("my_string","my string with spaces and @&$@(#!@&!"),
                                              ("a_field_size",2),
                                              ("float_array",(3.0,-4586.2)),
                                              ("alpha",set(list("a")))])),
                                 (("200    5, 212341,12341323    my string with spaces and @&$@(#!@&!    2    3.0,-4586.2    ".replace("    ","\t")),
                                 OrderedDict([("number",200),
                                              ("points",(5,212341,12341323)),
                                              ("my_string","my string with spaces and @&$@(#!@&!"),
                                              ("a_field_size",2),
                                              ("float_array",(3.0,-4586.2)),
                                              ("alpha",set())])),
                                ]]*3
