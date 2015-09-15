#!/usr/bin/env python
"""This module contains functions for escaping, unescaping, and parsing 
tokens from the ninth column of `GTF2`_ and `GFF3`_ files.

Important methods
-----------------
:py:func:`make_GTF2_tokens`
    Format a dictionary of attributes as `GTF2`_ column 9 attributes

:py:func:`make_GFF3_tokens`
    Format a dictionary of attributes as `GFF3`_ column 9 attributes

:py:func:`parse_GTF2_tokens`
    Parse `GTF2`_ column 9 tokens into a dictionary of key-value pairs

:py:func:`parse_GFF3_tokens`
    Parse `GFF3`_ column 9 tokens into a dictionary of key-value pairs

See also
--------
  - `The Sequence Ontology GFF3 specification <http://www.sequenceontology.org/gff3.shtml>`_
  - `The Brent lab GTF2.2 specification <http://mblab.wustl.edu/GTF22.html>`_
"""

# Unit tests for these are in :py:mod:`plastid.test.unit.genomics.readers.test_gff`



import re
import shlex
import copy
import warnings
from plastid.util.services.exceptions import FileFormatWarning

gtfpat = re.compile(r"^ *([^ ]*) +(.*) *$")

# From the spec: http://www.sequenceontology.org/gff3.shtml
# In addition to Parent, the Alias, Note, Dbxref and Ontology_term attributes can have multiple values.
# Also, SGD uses 'dbxref' instead of 'Dbxref'
_GFF3_DEFAULT_LISTS=("Parent","Alias","Note","Dbxref","Ontology_term","dbxref")


#===============================================================================
# INDEX: helper functions for escaping
#===============================================================================

# must escape % first, otherwise we'll end up escaping everything else, 
# since other escape codes start with percent signs
_GFF3_escape_sequences = [
 ('%', '%25'), # percent signs MUST be escaped FIRST 
 (';', '%3B'),
 (',', '%2C'),
 ('=', '%3D'),
 ('&', '%26'),
 ('\x00', '%00'),
 ('\x01', '%01'),
 ('\x02', '%02'),
 ('\x03', '%03'),
 ('\x04', '%04'),
 ('\x05', '%05'),
 ('\x06', '%06'),
 ('\x07', '%07'),
 ('\x08', '%08'),
 ('\t', '%09'),
 ('\n', '%0A'),
 ('\x0b', '%0B'),
 ('\x0c', '%0C'),
 ('\r', '%0D'),
 ('\x0e', '%0E'),
 ('\x0f', '%0F'),
 ('\x10', '%10'),
 ('\x11', '%11'),
 ('\x12', '%12'),
 ('\x13', '%13'),
 ('\x14', '%14'),
 ('\x15', '%15'),
 ('\x16', '%16'),
 ('\x17', '%17'),
 ('\x18', '%18'),
 ('\x19', '%19'),
 ('\x1a', '%1A'),
 ('\x1b', '%1B'),
 ('\x1c', '%1C'),
 ('\x1d', '%1D'),
 ('\x1e', '%1E'),
 ('\x1f', '%1F'),
 ('\x7f', '%7F'),
 ('\x80', '%80'),
 ('\x81', '%81'),
 ('\x82', '%82'),
 ('\x83', '%83'),
 ('\x84', '%84'),
 ('\x85', '%85'),
 ('\x86', '%86'),
 ('\x87', '%87'),
 ('\x88', '%88'),
 ('\x89', '%89'),
 ('\x8a', '%8A'),
 ('\x8b', '%8B'),
 ('\x8c', '%8C'),
 ('\x8d', '%8D'),
 ('\x8e', '%8E'),
 ('\x8f', '%8F'),
 ('\x90', '%90'),
 ('\x91', '%91'),
 ('\x92', '%92'),
 ('\x93', '%93'),
 ('\x94', '%94'),
 ('\x95', '%95'),
 ('\x96', '%96'),
 ('\x97', '%97'),
 ('\x98', '%98'),
 ('\x99', '%99'),
 ('\x9a', '%9A'),
 ('\x9b', '%9B'),
 ('\x9c', '%9C'),
 ('\x9d', '%9D'),
 ('\x9e', '%9E'),
 ('\x9f', '%9F')
 ]
"""List mapping characters to their escape sequences, per the `GFF3`_ specification"""

_GTF2_escape_sequences = copy.deepcopy(_GFF3_escape_sequences)
_GTF2_escape_sequences.append(("\"","%22"))
"""List mapping characters to their escape sequences for `GTF2`_. These are undefined,
but we are using `GFF3`_ characters plus double quotation marks as a convention.
"""


def escape(inp,char_pairs):
    """Escape reserved characters specified in the list of tuples `char_pairs`
    
    Parameters
    ----------
    inp : str
        Input string
    
    chair_pairs : list
        List of tuples of (character, escape sequence for character)
    
    
    Returns
    -------
    str
        Escaped output
    
    
    See also
    --------
    unescape_GFF3
    """
    for char_, repl in char_pairs:
        inp = inp.replace(char_,repl)
    
    return inp

def unescape(inp,char_pairs):
    """Unescape reserved characters specified in the list of tuples `char_pairs`
    
    Parameters
    ----------
    inp : str
        Input string
    
    
    Returns
    -------
    str
        Unescaped output
    
    
    See also
    --------
    escape_GFF3
    """   
    for repl, char_ in reversed(char_pairs): # reverse order of escaping/unescaping
        inp = inp.replace(char_,repl)
    
    return inp

def escape_GFF3(inp):
    """Escape reserved characters in `GFF3`_ tokens using percentage notation.
    
    In the `GFF3`_ spec, reserved characters include:
    
        - control characters (ASCII 0-32, 127, and 128-159)
        
        - tab, newline, & carriage return
        
        - semicolons & commas
        
        - the percent sign
        
        - the equals sign
        
        - the ampersand  
    
    Parameters
    ----------
    inp : str
        Input string
    
    chair_pairs : list
        List of tuples of (character, escape sequence for character)
    
    
    Returns
    -------
    str
        Escaped output
    
    
    See also
    --------
    unescape_GFF3
    """    
    return escape(inp,_GFF3_escape_sequences)

def unescape_GFF3(inp):
    """Unescape reserved characters in `GFF3`_ tokens using percentage notation.
    
    In the `GFF3`_ spec, reserved characters include:
    
        - control characters (ASCII 0-32, 127, and 128-159)
        
        - tab, newline, & carriage return
        
        - semicolons & commas
        
        - the percent sign
        
        - the equals sign
        
        - the ampersand  
    
    Parameters
    ----------
    inp : str
        Input string
    
    
    Returns
    -------
    str
        Unescaped output
    
    
    See also
    --------
    escape_GFF3
    """    
    return unescape(inp,_GFF3_escape_sequences)

def escape_GTF2(inp):
    """Escape reserved characters in `GTF2`_ tokens using percentage notation.
    While the `GTF2`_ spec is agnostic for escaping, it is useful when adding
    extra attributes to files. As a convention, we escape the characters
    specified in the `GFF3`_ spec, as well as double quotation marks.
    
    In the `GTF2`_ spec, reserved characters include:
    
        - control characters (ASCII 0-32, 127, and 128-159)
        
        - tab, newline, & carriage return
        
        - semicolons & commas
        
        - the percent sign
        
        - the equals sign
        
        - the ampersand  
    
    Parameters
    ----------
    inp : str
        Input string
    
    chair_pairs : list
        List of tuples of (character, escape sequence for character)
    
    
    Returns
    -------
    str
        Escaped output
    
    
    See also
    --------
    unescape_GFF3
    """    
    return escape(inp,_GTF2_escape_sequences)

def unescape_GTF2(inp):
    """Unescape reserved characters in `GTF2`_ tokens using percentage notation.
    While the `GTF2`_ spec is agnostic for escaping, it is useful when adding
    extra attributes to files. As a convention, we escape the characters
    specified in the `GFF3`_ spec, as well as single quotation marks.
        
    In the `GFF3`_ spec, reserved characters include:
    
        - control characters (ASCII 0-32, 127, and 128-159)
        
        - tab, newline, & carriage return
        
        - semicolons & commas
        
        - the percent sign
        
        - the equals sign
        
        - the ampersand  
    
    Parameters
    ----------
    inp : str
        Input string
    
    
    Returns
    -------
    str
        Unescaped output
    
    
    See also
    --------
    escape_GFF3
    """    
    return unescape(inp,_GTF2_escape_sequences)

#===============================================================================
# INDEX: attribute token formatting and parsing
#===============================================================================

def _make_generic_tokens(attr,excludes=[],join_pat='%s %s; ',escape=None):
    """Helper function to convert the `attr` dict of a |SegmentChain|
    into the string representation used in GFF files. This includes
    URL escaping of keys and values, and catenating lists with `','`
    before string conversion
    
    Parameters
    ----------
    attr : dict
        Dictionary of key-value pairs to export
        
    excludes : list<str>
        List of keys to exclude from string
        
    join_pat
        printf-style pattern explaining how to join key:value pairs
    
    escape : None or func, optional
        If None, no special characters are escaped. If a function, that
        funciton will be used to perform the escaping. (Default: `False`)
        
    Returns
    -------
    str
    """
    f = lambda x: x[0] not in excludes
    if escape is None:
        esc = lambda inp: inp
    else:
        esc = lambda inp: escape(str(inp))
    ltmp = []
    for key, val in filter(f,attr.items()):
        if isinstance(val,list):
            val = ",".join([esc(X) for X in val])
        else:
            val = esc(val)
        ltmp.append(join_pat % (esc(key),val))
    #ltmp = [join_pat % (esc(key),esc(val)) for (key,val) in filter(f,attr.items())]
    return ''.join(ltmp)

def make_GFF3_tokens(attr,excludes=[],escape=True):
    """Helper function to convert the `attr` dict of a |SegmentChain|
    into the string representation used in `GFF3`_ files. This includes
    URL escaping of special characters, and catenating lists with '`,`'
    before string conversion

    Examples
    --------
        >>> d = {'a':1,'b':2,'c':3,'d':4,'e':5,'z':26,'text':"something; with escape sequences"}
        >>> _make_GFF3_tokens(d)
        'a=1;c=3;b=2;e=5;d=4;z=26;text=something%3B with escape sequences'
    
        >>> excludes=['a','b','c']
        >>> _make_GFF3_tokens(d,excludes)
        'e=5;d=4;z=26;text=something%3B with escape sequences'

        >>> d = {'a':1,'b':2,'c':[3,7],'d':4,'e':5,'z':26}
        >>> _make_GFF3_tokens(d)
        'a=1;c=3,7;b=2;e=5;d=4;z=26'


    Parameters
    ----------
    attr : dict
        Dictionary of key-value pairs to export
        
    excludes : list<str>
        List of keys to exclude from string
        
    escape : bool
        If True, special characters in output are `GFF3`_-escaped (Default: `True`)
        
    Returns
    -------
    str
        Data formatted for *attributes* column of `GFF3`_ (column 9)
    """
    if escape == True:
        escape = escape_GFF3
    else:
        escape = None
    
    return _make_generic_tokens(attr,excludes=excludes,join_pat="%s=%s;",escape=escape)

def make_GTF2_tokens(attr,excludes=[],escape=True):
    """Helper function to convert the `attr` dict  of a |SegmentChain|
    into the string representation used in `GTF2`_ files. By default, special
    characters defined in the `GFF3`_ spec will be URL-escaped.

    Examples
    --------
        >>> d = {'transcript_id' : 't;id', 'a':1,'b':2,'c':3,'d':4,'e':5,'z':26,
                    'gene_id' : 'gid'}
        >>> _make_GTF2_tokens(d)
        'transcript_id "t%3Bid"; gene_id "gid"; a "1"; c "3"; b "2"; e "5"; d "4"; z "26";'
    
        >>> excludes=['a','b','c']
        >>> _make_GTF2_tokens(d,excludes)
        'transcript_id "t%3Bid"; gene_id "gid"; e "5"; d "4"; z "26";'


    Parameters
    ----------
    attr : dict
        Dictionary of key-value pairs to export
        
    excludes : list<str>
        List of keys to exclude from string
        
    escape : bool
        If True, special characters in output are `GTF2`_-escaped (Default: `True`)
        
    Returns
    -------
    str
        Data formatted for *attributes* column of `GTF2`_ (column 9)
    """
    excludes.extend(["transcript_id","gene_id"])    
    stmp = 'gene_id "%s"; transcript_id "%s"; ' % (attr.get("gene_id"),
                                                   attr.get("transcript_id")) 
    
    if escape == True:
        escape = escape_GTF2
    else:
        escape = None

    return stmp + _make_generic_tokens(attr,excludes=excludes,join_pat='%s "%s"; ',escape=escape).strip(" ")

def parse_GFF3_tokens(inp,list_types=_GFF3_DEFAULT_LISTS):
    """Helper function to parse tokens in the final column of a `GFF3`_ file
    into a dictionary of attributes. Because, the following attributes are
    permitted to have multiple values in the `GFF3`_ spec, their values, if present
    are returned as lists in the dictionary rather than strings:
    
        - `Parent`
        - `Alias`
        - `Note`
        - `Dbxref`
        - `Ontology_term`
 
    All values are unescaped folowing the `GFF3`_ specification.
 
    Examples
    --------
        >>> tokens = 'a=1;c=3;b=2;e=5;d=4;z=26,Parents=gene01'
        >>> parse_GFF3_tokens(tokens)
        {'a': '1', 'c': '3', 'b': '2', 'e': '5', 'd': '4', 'z': '26', 'parents' : ['gene01'] }

        >>> tokens = 'a=1;c=3,7;b=2;e=5;d=4;z=26,Parents=gene01,gene02'
        >>> parse_GFF3_tokens(tokens)
        {'a': '1', 'c': '3,7', 'b': '2', 'e': '5', 'd': '4', 'z': '26', 'parents' : ['gene01','gene02']}

 
    Parameters
    ----------
    inp : str
        Ninth column of `GFF3`_ entry
    
    list_types : str
        Names of attributes that should be returned as lists
        (Default: %s)
         
    Returns
    -------
    dict : key-value pairs
    """ % ",".join(_GFF3_DEFAULT_LISTS)
    d = {}
    items = inp.strip("\n").strip(";").split(";")
    for item in items:
        if len(item) > 0:
            key, val = item.split("=")
            key = unescape_GFF3(key.strip(" "))
            if key in list_types:
                val = [unescape_GFF3(X) for X in val.strip(" ").split(",")]
            else:
                val = unescape_GFF3(val.strip(" "))
                
            if key in d:
                warnings.warn("Found duplicate attribute key '%s' in GFF3 line. Catenating value with previous value for key in attr dict:\n    %s" % (key,inp),
                              FileFormatWarning)
                val = "%s,%s" % (d[key],val)
            d[key] = val
    return d

def parse_GTF2_tokens(inp):
    """Helper function to parse tokens in the final column of a `GTF2`_ file
    into a dictionary of attributes. All attributes are returned as strings,
    and are unescaped if GFF escape sequences (e.g. *'%2B'*) are present.

    If duplicate keys are present (e.g. as in GENCODE `GTF2`_ files),
    their values are catenated, separated by a comma.
    
    Examples
    --------
        >>> tokens = 'gene_id "mygene"; transcript_id "mytranscript";'
        >>> parse_GTF2_tokens(tokens)
        {'gene_id' : 'mygene', 'transcript_id' : 'mytranscript'}
    
        >>> tokens = 'gene_id "mygene"; transcript_id "mytranscript"'
        >>> parse_GTF2_tokens(tokens)
        {'gene_id' : 'mygene', 'transcript_id' : 'mytranscript'}
    
        >>> tokens = 'gene_id "mygene;"; transcript_id "myt;ranscript"'
        >>> parse_GTF2_tokens(tokens)
        {'gene_id' : 'mygene;', 'transcript_id' : 'myt;ranscript'}
    
        >>> tokens = 'gene_id "mygene"; transcript_id "mytranscript"; tag "tag value";'
        >>> parse_GTF2_tokens(tokens)
        {'gene_id' : 'mygene', 'tag' : 'tag value', 'transcript_id' : 'mytranscript'}

        >>> tokens = 'gene_id "mygene"; transcript_id "mytranscript"; tag "tag value"; tag "tag value 2";'
        >>> parse_GTF2_tokens(tokens)
        {'gene_id' : 'mygene', 'tag' : 'tag value,tag value 2', 'transcript_id' : 'mytranscript'}



    Parameters
    ----------
    inp : str
        Ninth column of `GTF2`_ entry
        
    Returns
    -------
    dict : key-value pairs
    """
    d = {}
    items = shlex.split(inp.strip("\n"))
    assert len(items) % 2 == 0
    for i in range(0,len(items),2):
        key = unescape_GTF2(items[i])
        val = items[i+1]
        # require separation by semicolons for all but final token
        if i+1 < len(items) - 2:
            assert val.endswith(";")
        
        if val.endswith(";"):
            val = val[:-1]

        if key in d:
            warnings.warn("Found duplicate attribute key '%s' in GTF2 line. Catenating value with previous value for key in attr dict:\n    %s" % (key,inp),
                          FileFormatWarning)
            d[key] = "%s,%s" % (d[key],unescape_GTF2(val))

        else:
            d[key] = unescape_GTF2(val)
        
    return d
