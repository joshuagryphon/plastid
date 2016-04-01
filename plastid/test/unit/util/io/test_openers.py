#!/usr/bin/env python
from plastid.util.io.openers import (get_short_name, pretty_print_dict)
from nose.tools import assert_equal

def test_get_short_name():
    tests = [("test","test",{}),
             ("test.py","test",dict(terminator=".py")),
             ("/home/jdoe/test.py","test",dict(terminator=".py")),
             ("/home/jdoe/test.py.py","test.py",dict(terminator=".py")),
             ("/home/jdoe/test.py.2","test.py.2",{}),
             ("/home/jdoe/test.py.2","test.py.2",dict(terminator=".py")),
             ("plastid.bin.test","test",dict(separator="\.",terminator=""))
             ]
    for inp, expected, kwargs in tests:
        found = get_short_name(inp,**kwargs)
        msg = "test_get_short_name(): failed on input '%s'. Expected '%s'. Got '%s'" % (inp,expected,found)
        yield assert_equal, expected, found, msg

    
def test_pretty_print_dict():
    dtmp = { "a" : 1,
             "b" : 2.3,
             "c" : "some string",
             "d" : "some string with 'subquotes' inside",
             "e" : (3,4,5),
             "somereallyreallylongname" : "short val",
            }
    expected = """{
          'a'                        : 1,
          'b'                        : 2.3,
          'c'                        : 'some string',
          'd'                        : 'some string with 'subquotes' inside',
          'e'                        : (3, 4, 5),
          'somereallyreallylongname' : 'short val',
}
"""
    found = pretty_print_dict(dtmp)
    assert_equal(expected,found,"Dictionary did not pretty-print!\nExpected:\n%s\n\nFound:\n%s\n\n" % (expected,found)) 