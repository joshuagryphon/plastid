#!/usr/bin/env python
"""Tools for converting declarations written in Jim Kent's & Heidi Brumbaugh's
`autoSql <http://www.linuxjournal.com/article/5949>`_ table/object specification language
into Python objects that can parse blocks of delimited text into records that follow
the specification of the autoSql declaration, each field in each record converted
to its native type.

The classes here are initialized with autoSql declarations. Then, when instances
are called on blocks of text, they parse those blocks according to the specification
of the declaration.

Examples
--------
Create a parser for some data structure based upon an autoSql declaration::

    >>> declaration = '''table easy_table
    "A table with a comment on the next line" 
        (
        uint number auto; "a number with a token"
        uint [3] points ; "r,g,b values"
        lstring  my_string ; "a long string"
        uint a_field_size ; "the size for the next field"
        float [a_field_size] float_array ; "an array of floats"
        set(a,b,c) alpha ; "the first three letters of the alphabet"
        )
    '''
    >>> record_parser = AutoSqlDeclaration(declaration)


Parse text into a record using that parser::

    >>> record_parser.parse("3    1,2,3    my string with spaces    5    1.1,1.2,1.3,1.4,1.5    a,b")
    OrderedDict([("number",3),
                 ("points",(1,2,3)),
                 ("my_string","my string with spaces"),
                 ("a_field_size",5),
                 ("float_array",(1.1,1.2,1.3,1.4,1.5)),
                 ("alpha",{'a','b'}]))


Important classes
-----------------
|AutoSqlDeclaration|
    Parses autoSql declarations for *table, simple,* and *object* types,
    and delegates parsing of its own fields to appropriate parsers 
    (e.g. |AutoSqlField|, |SizedAutoSqlField|, and others).

|AutoSqlField|, |SizedAutoSqlField|, |ValuesAutoSqlField|, et c
    Parse various sorts of fields within an autoSql declaration block


Notes
-----
1. These parsers seek only to provide Python bindings for autoSql declarations,
**NOT** to convert autoSql to C or SQL, functions which are already provided
by `Jim Kent's utilities <https://github.com/ENCODE-DCC/kentUtils/tree/master/>`_

2. ``set`` and ``enum`` field types are both returned as ``set`` s of strings

3. ``primary``, ``index``, and ``auto`` SQL tags are accepted in input, but ignored,
because they are not relevant for our purposes

4. It is assumed, as would be the case for BigBed files, that text blocks passed
to the ``__call__()`` methods of |AutoSqlDeclaration| objects
are tab-delimited. This behavior can be controlled by setting the ``delim``
attribute of the various parsers

5. Although declarations are routinely nested as fields within other
declarations in C ``struct`` s and in SQL databases, in the absence of a standard,
it is unclear how these would be serialized within tab-delimited columns of `BigBed`_
or BigWig files. Therefore, we do not support this. If or when such a standard is
described, we will implement it.


See Also
--------
`Kent & Brumbaugh, 2002 <http://www.linuxjournal.com/article/5949>`_
    First publication of autoSql & autoXml 

`Updated autoSql Grammar specification <https://github.com/ENCODE-DCC/kentUtils/blob/36d6274459f644d5400843b8fa097b380b8f7867/src/hg/autoSql/autoSql.doc>`_
    Explanation of autoSql grammar

The ENCODE project's `tests for autoSql parsers <https://github.com/ENCODE-DCC/kentUtils/tree/master/src/hg/autoSql/tests/input>`_
    Official autoSql unit tests
"""
import re
import warnings
from collections import OrderedDict
from abc import abstractmethod

# regular expressions that recognize various autoSql elements
_pattern_bits = { "start"   : r"^\s*",
                  "type"    : r"(?P<type>\w+)",
                  "name"    : r"\s+(?P<name>\w+)\s*",
                  "semi"    : r"\s*;\s*",
                  "comment" : r"\"(?P<comment>[^\"]*)\"",
                  "size"    : r"\s*\[\s*(?P<size>\w+)\s*\]\s*",
                  "values"  : r"\s*\(\s*(?P<value_names>[^()]+)\s*\)\s*",
                  "optionals" : r"(?P<opt1>\s+primary|\s+auto|\s+index\s*(\[\s*\d+\s*\])?)?\s*(?P<opt2>\s+primary|\s+auto|\s+index\s*(\[\s*\d+\s*\])?)?\s*(?P<opt3>\s+primary|\s+auto|\s+index\s*(\[\s*\d+\s*\])?)?", 
                  "declare_type_name" : r"(?P<declare_type>object|simple|table)\s+(?P<declare_name>\w+)\s+",
                  #"field_text" :  r"\s*\(\s*(?P<field_text>.+)\s*\)",
                  "field_text" :  r"\s*\(\s*(?P<field_text>.*)\)",
                 }


class AbstractAutoSqlElement(object):
    """Abstract base class for parsers of autoSql elements
    
    Attributes
    ----------
    attr : dict
        Dictionary of attributes describing the element (e.g. *name,* *type,* et c) 
        
    autosql : str
        Block of autoSql text specifying format of element
        
    match_pattern : :py:class:`re.RegexObject`
        Pattern that determines whether or not a block of autoSql matches this object
    
    parent : instance of subclass of |AbstractAutoSqlElement|, or None
        Parent / enclosing element
    
    field_types : dict
        Dictionary matching type names (as strings) to formatters that parse them
        from plaintext
    
    delim : str, optional
        Text delimiter for fields in blocks called by :py:meth:~__call__~
        (Default: "\t")
    """
    match_str = ""
    match_pattern = re.compile(match_str)
     
    def __init__(self,autosql,parent=None,delim="\t"):
        self.autosql = autosql
        self.parent  = parent
        self.delim   = delim
        self.field_types = { "int"    : (int,    "i"),  #32-bit
                             "uint"   : (int,    "I"), #32-bit
                             "short"  : (int,    "h"), #16-bit
                             "ushort" : (int,    "H"), #16-bit
                             "byte"   : (int,    "b"), #8-bit
                             "ubyte"  : (int,    "B"), #8-bit
                             "float"  : (float,  "f"), #single-precision
                             "char"   : (str,    "c"), #8-bit
                             "string" : (str,    "s"), #variable up to 255bytes
                             "lstring": (str,    "s"),  #variable up to 2billion bytes
                           }
        self.attr = self.match_pattern.search(autosql).groupdict()
    
    def __repr__(self):
        return "<%s name=%s type=%s>" % (self.__class__.__name__,
                                         self.attr["name"],
                                         self.attr.get("type",self.__class__.__name__))
    
    def add_type(self,name,formatter):
        """Add a type to the parser
        
        Parameters
        ----------
        name : str
            Name of data type
        
        formatter : callable
            Function/callable that, when applied to autoSql text, yields
            an object of the type specified by ``name``
        """
        self.field_types[name] = formatter
        
    @abstractmethod
    def __call__(self,text,rec=None):
        """Parse an OrderedDict matching ``self.autosql`` from a block of delimited text
        
        Parameters
        ----------
        text : str
            Multiline text block, formatted in autoSql
        
        rec : OrderedDict or None, optional
            Record whose attributes are being populated by recursive
            processing of ``text``. Passed in cases where fields sized by variables
            need to look up instance values of earlier fields to evaluate those
            variables.
        """
        pass
    
    @staticmethod
    def mask_comments(text):
        """Mask all comments in an autoSql block in order to facilitate parsing
        by regular expressions
        
        Parameters
        ----------
        text : str
            autoSql-formatted text
        
        Returns
        -------
        str
            Text with comments replaced by "xxxxxx" of same length
        
        list
            List of (comment.start,comment.end), including quotes, for each comment
            in ``text`` 
        """
        cpat = re.compile(r"\"[^\"]+\"")
        match_locs = []
        for match in cpat.finditer(text):
            my_start = match.start()
            my_end   = match.end()
            match_len = my_end - my_start
            match_locs.append((my_start+1,my_end-1))
            text = text[:my_start+1] + "x"*(match_len-2) + text[my_end-1:]
        
        return text, match_locs

    @classmethod
    def matches(cls,text):
        """Determine whether autoSql formatting text matches this autoSql element
        
        Parameters
        ----------
        text : str
            Block of autoSql-formatted declaration text
        
        Returns
        bool
            True an autoSql parser of this class's type can be made from this
            specification, otherwise False
        """
        return cls.match_pattern.search(text) is not None


class AutoSqlDeclaration(AbstractAutoSqlElement):
    """Parser factory that converts delimited text blocks into OrderedDicts,
    following the field names and types described by an autoSql declaration element
    
    Attributes
    ----------
    attr : dict
        Dictionary of descriptive attributes (e.g. *name,* *type,* *declare_type,* et c) 
    
    field_formatters : OrderedDict
        Dictionary mapping field names to type names

    field_comments : OrderedDict
        Dictionary mapping field names to comments
        
    field_types : dict
        Dictionary matching type names (as strings) to formatters that parse them
        from plaintext
    
    autosql : str
        Block of autoSql text specifying format of element
        
    match_pattern : :py:class:`re.RegexObject`
        Pattern that determines whether or not a block of autoSql matches this object
    
    parent : instance of subclass of |AbstractAutoSqlObject|, or None
        Parent / enclosing element. Default: None
    
    delim : str, optional
        Text delimiter for fields in blocks called by :py:meth:~__call__~
        (Default: "\t")
    
    
    Methods
    -------
    :py:meth:`AutoSqlDeclaration.__call__`
        Parse autoSql-formatted blocks of text according to this declaration
    """
    
    match_str  = r"".join([_pattern_bits[X] for X in ("start","declare_type_name","comment","field_text")])
    match_pattern = re.compile(match_str,re.S)

    def __init__(self,autosql,parent=None,delim="\n"):
        """Create an |AutoSqlDeclaration|
        
        Parameters
        ----------
        autosql : str
            Block of autoSql text specifying format of element
            
        parent : instance of subclass of |AbstractAutoSqlObject| or None, optional
            Parent / enclosing element. Default: None
        
        delim : str, optional
            Field delimiter (default: tab)
        """
        AbstractAutoSqlElement.__init__(self,autosql,parent=parent,delim="\t")
        
        # re-do regex match masking out comments, in case the comments
        # contain special characters that would mess up the parsing
        masked_sql, comment_match_locs = self.mask_comments(autosql)
        match_dict = self.match_pattern.search(masked_sql).groupdict()

        self.attr["declare_type"] = match_dict["declare_type"]
        self.attr["name"]         = match_dict["declare_name"]
        masked_field_text = match_dict["field_text"]

        self.attr["comment"] = autosql[comment_match_locs[0][0]:comment_match_locs[0][1]].strip("\n").strip("\"")
        field_text_start = masked_sql.index(masked_field_text)
        self._field_text  = autosql[field_text_start:field_text_start+len(masked_field_text)]
        
        if self.parent is not None:
            self.parent.add_type(self.attr["declare_name"],self)
            
        self.field_formatters = OrderedDict()
        self.field_comments   = OrderedDict()
        self._parse_fields()

    def _parse_fields(self):
        """Parse fields of an autoSql declaration, and populate
        ``self.field_formatters`` and ``self.field_comments``.
        """
        # order in which we try to match autoSql fields        
        match_order = [AutoSqlField,SizedAutoSqlField,ValuesAutoSqlField]

        # fields are area of string from last starting point to end of comment
        # first starting point is 0;all subsequent starting points will be end 
        # of previous comment
        
        _, comment_locs = self.mask_comments(self._field_text)
        last_index = 0
        for (_,next_index) in comment_locs:
            field_str = self._field_text[last_index:next_index+1]
            for field_class in match_order:
                if field_class.matches(field_str):
                    my_parser = field_class(field_str)
                    self.field_formatters[my_parser.attr["name"]]  = my_parser
                    self.field_comments[my_parser.attr["name"]]    = my_parser.attr["comment"]
            
            last_index = next_index+1

    def __repr__(self):
        return "<%s name=%s type=%s fields=[%s]>" % (self.__class__.__name__,
                                         self.attr["name"],
                                         self.attr.get("type",self.__class__.__name__),
                                         ",".join(self.field_formatters.keys()))
        
    def __call__(self,text,rec=OrderedDict()):
        """Parse an OrderedDict matching ``self.autosql`` from a block of delimited text
        
        Parameters
        ----------
        text : str
            Multiline text block, formatted in autoSql

        rec : OrderedDict or None, optional
            Record whose attributes are being populated by recursive
            processing of ``text``. Passed in cases where fields sized by variables
            need to look up instance values of earlier fields to evaluate those
            variables.
        
        Returns
        -------
        OrderedDict
            Dictionary mapping field names to their values
        """
        items = text.split(self.delim)
        obj = OrderedDict()
        for item, (field_name,formatter) in zip(items,self.field_formatters.items()):
            obj[field_name] = formatter(item,rec=obj)
        
        return obj


class AutoSqlField(AbstractAutoSqlElement):
    """Parser factory for autoSql fields of type ``fieldType fieldName ';' comment``
    
    Attributes
    ----------
    attr : dict
        Dictionary of descriptive attributes (e.g. *name,* *type,* et c) 

    formatter : callable
        Callable/function that converts plain text into an object of the correct type
        
    autosql : str
        Block of autoSql text specifying format of element
        
    match_pattern : :py:class:`re.RegexObject`
        Pattern that determines whether or not a block of autoSql matches this object
    
    parent : instance of subclass of |AbstractAutoSqlObject|, or None
        Parent / enclosing element. Default: None
    
    delim : str, optional
        Text delimiter for fields in blocks called by :py:meth:~__call__~
        (Default: "\n")
        
    Methods
    -------
    :py:meth:`AutoSqlField.__call__`
        Parse autoSql-formatted blocks of text into the object type specified
        by this field        
    """
    match_str = r"".join([_pattern_bits[X] for X in ("start","type","name","optionals","semi","comment")])
    match_pattern = re.compile(match_str)

    def __init__(self,autosql,parent=None,delim=""):
        """Create an |AutoSqlField|
        
        Parameters
        ----------
        autosql : str
            Block of autoSql text specifying format of element
            
        parent : instance of subclass of |AbstractAutoSqlObject| or None, optional
            Parent / enclosing element. Default: None
        
        delim : str, optional
            Field delimiter (default: tab)
        """        
        AbstractAutoSqlElement.__init__(self,autosql,parent=parent,delim=delim)
        try:
            self.formatter = self.field_types[self.attr["type"]][0]
        except KeyError:
            self.formatter = self.parent.field_types[self.attr["type"]][0]
    
    def __call__(self,text,rec=None):
        """Parse an value matching the field described by ``self.autosql``
        from a block of delimited text
        
        Parameters
        ----------
        text : str
            Multiline text block, formatted in autoSql
        
        Returns
        -------
        Value or object of appropriate type
        """
        try:
            return self.formatter(text)
        except ValueError:
            message = "Could not convert autoSql value '%s' for field '%s' to type '%s'. Leaving as str " % (text,
                                                                                                             self.attr["name"],
                                                                                                             self.formatter.__name__)
            warnings.warn(message,UserWarning) 
            return text


class SizedAutoSqlField(AutoSqlField):
    """Parser factory for autoSql fields of type ``fieldType `[` fieldSize `]` fieldName ';' comment``
    
    Attributes
    ----------
    attr : dict
        Dictionary of descriptive attributes (e.g. *name,* *size,* *type,* et c) 
    
    formatter : callable
        Callable/function that converts plain text into an object of the correct type
    
    autosql : str
        Block of autoSql text specifying format of element
        
    match_pattern : :py:class:`re.RegexObject`
        Pattern that determines whether or not a block of autoSql matches this object
    
    parent : instance of subclass of |AbstractAutoSqlObject|, or None
        Parent or enclosing element. Default: None
    
    delim : str, optional
        Text delimiter for fields in blocks called by :py:meth:~__call__~
        (Default: "\n")

    Methods
    -------
    :py:meth:`SizedAutoSqlField.__call__`
        Parse autoSql-formatted blocks of text into the tuples of the object type
        specified by this field  
    """    
    match_str = r"".join([_pattern_bits[X] for X in ("start","type","size","name","optionals","semi","comment")])
    match_pattern = re.compile(match_str)

    def __init__(self,autosql,size=1,parent=None,delim=","):
        """Create a |SizedAutoSqlField|
        
        Parameters
        ----------
        autosql : str
            Block of autoSql text specifying format of element
            
        parent : instance of subclass of |AbstractAutoSqlObject| or None, optional
            Parent / enclosing element. Default: None
        
        delim : str, optional
            Field delimiter (default: tab)
        """           
        AutoSqlField.__init__(self,autosql,parent=parent,delim=delim)
        try:
            self.attr["size"] = int(self.attr["size"])
            self.attr["size_is_int"] = True
        except ValueError:
            self.attr["size_is_int"] = False
    
    def __call__(self,text,rec=None):
        """Parse an value matching the field described by ``self.autosql``
        from a block of delimited text
        
        Parameters
        ----------
        text : str
            Multiline text block, formatted in autoSql

        rec : OrderedDict or None, optional
            Record whose attributes are being populated by recursive
            processing of ``text``. Passed in cases where fields sized by variables
            need to look up instance values of earlier fields to evaluate those
            variables.
        
        Returns
        -------
        tuple
            Tuple of appropriate type
        """
        if self.formatter != str:
            try:
                retval = tuple([self.formatter(X) for X in text.strip().strip(self.delim).split(self.delim)])
            except ValueError:
                message = "Could not convert autoSql value '%s' in field '%s' to tuple of type '%s'. Leaving as str " % (text,
                                                                                                                         self.attr["name"],
                                                                                                                         self.formatter.__name__)
                warnings.warn(message,UserWarning) 
                return text
        else:
            retval = text
        
        if self.attr["size_is_int"] == True:    
            assert len(retval) == self.attr["size"]
        else:
            assert len(retval) == rec[self.attr["size"]]
        
        return retval


# for set, enum types
class ValuesAutoSqlField(AbstractAutoSqlElement):
    """Parser factory for autoSql fields of type ``fieldType `(` fieldValues `)` fieldName ';' comment``
    where ``fieldType`` would typically be ``set`` or ``enum``
    """
    
    match_str = r"".join([_pattern_bits[X] for X in ("start","type","values","name","optionals","semi","comment")])
    match_pattern = re.compile(match_str)
    
    def __init__(self,autosql,parent=None,delim=","):
        """Create a |ValuesAutoSqlField|
        
        Parameters
        ----------
        autosql : str
            Block of autoSql text specifying format of element
            
        parent : instance of subclass of |AbstractAutoSqlObject| or None, optional
            Parent / enclosing element. Default: None
        
        delim : str, optional
            Field delimiter (default: tab)
        """            
        AbstractAutoSqlElement.__init__(self,autosql,parent=parent,delim=delim)
        self.attr["value_names"] = [X.strip() for X in self.attr["value_names"].split(",")]

    def __call__(self,text,rec=None):
        """Parse an value matching the field described by ``self.autosql``
        from a block of delimited text
        
        Parameters
        ----------
        text : str
            Multiline text block, formatted in autoSql

        rec : OrderedDict or None, optional
            Record whose attributes are being populated by recursive
            processing of ``text``. Passed in cases where fields sized by variables
            need to look up instance values of earlier fields to evaluate those
            variables.
        
        Returns
        -------
        set
            set of items found in column 
        """
        items = set([X.strip() for X in text.strip(self.delim).split(self.delim) if len(X.strip()) > 0])
        return items


# @notimplemented
# class DeclareTypeAutoSqlField(AbstractAutoSqlElement):
#     """Parse an autoSql field of type ``declareType declareName fieldName ';' comment``"""
#     match_str = r"".join([_pattern_bits[X] for X in ("declare_type_name","name","semi","comment")])
#     match_pattern = re.compile(match_str)
#     
#     def __init__(self,autosql,parent=None,delim=" "):
#         AbstractAutoSqlElement.__init__(self,autosql,parent=parent,delim=delim)
#         
#         if self.parent is not None:
#             self.parent.add_type(self.attr["declare_name"],self)
# 
#     def __call__(self,text,rec=None):
#         """Parse this type of object from a block of autoSql formatted text
#         
#         Parameters
#         ----------
#         text : str
#             Multiline text block, formatted in autoSql
#         
#         rec : OrderedDict or None, optional
#             Record whose attributes are being populated by recursive
#             processing of ``text``. Passed in cases where fields sized by variables
#             need to look up instance values of earlier fields to evaluate those
#             variables.
#         """
#         pass
#         
#         
# @notimplemented
# class SizedDeclareTypeAutoSqlField(DeclareTypeAutoSqlField):
#     """Parse an autoSql field of type ``declareType declareName `[` fieldSize `]` fieldName ';' comment``"""
#     
#     match_str = r"".join([_pattern_bits[X] for X in ("declare_type_name","size","name","semi","comment")])
#     match_pattern = re.compile(match_str)
# 
#     def __call__(self,text,rec=None):
#         """Parse this type of object from a block of autoSql formatted text
#         
#         Parameters
#         ----------
#         text : str
#             Multiline text block, formatted in autoSql
#         
#         rec : OrderedDict or None, optional
#             Record whose attributes are being populated by recursive
#             processing of ``text``. Passed in cases where fields sized by variables
#             need to look up instance values of earlier fields to evaluate those
#             variables.
#         """
#         pass
# 
# 
# @notimplemented
# class ValuesDeclareTypeAutoSqlField(DeclareTypeAutoSqlField):
#     """Parse an autoSql field of type ``declareType declareName `(` fieldValues `)` fieldName ';' comment``"""
#     
#     match_str = r"".join([_pattern_bits[X] for X in ("declare_type_name","values","name","semi","comment")])
#     match_pattern = re.compile(match_str)
#     
#     def __call__(self,text,rec=None):
#         """Parse this type of object from a block of autoSql formatted text
#         
#         Parameters
#         ----------
#         text : str
#             Multiline text block, formatted in autoSql
#         
#         rec : OrderedDict or None, optional
#             Record whose attributes are being populated by recursive
#             processing of ``text``. Passed in cases where fields sized by variables
#             need to look up instance values of earlier fields to evaluate those
#             variables.
#         """
#         pass   
