#!/usr/bin/env python
"""Tools for reading values from binary files

See Also
--------
:py:mod:`struct`
    Binary data structures in Python
"""
import numpy
import sys
import struct
from collections import namedtuple

class BinaryParserFactory(object):
    """Parser factory for different types of binary records.
    
    Creates parsers that unpack binary byte streams into dictionaries
    that match field names to values. These parsers are most useful as components
    of binary file readers.
    
    Attributes
    ----------
    name : str
        Human-readable name for parser
    
    fmt : str
        String specifying binary format of data, as specified in :py:mod:`struct`
    
    fields : list
        List of strings specifying variable names to bind to data
        when unpacked from a binary file, in same order as items in ``fmt``
    
    nt : :class:`~collections.namedtuple`
        A :class:`~collections.namedtuple` instance that will provide names
        to the unpacked data

    
    Examples
    --------
    A binary RGB color parser::
    
        >>> ColorParser = BinaryParserFactory("ColorParser","3Q",["r","g","b"])
        >>> fh = open("some_binary_file_containing_colors.bin","rb") # 'b' is important in mode flag!!
        >>> fh.seek(byte_location_of_an_rgb_color)
        >>> rgb_dict = ColorParser(fh) # read and parse 3 8-bit integers from file
        >>> rgb_dict
            { "r" : 255,
              "g" : 0,
              "b" : 52 }


    See Also
    --------
    struct
        For information on format strings
    """
    
    def __init__(self,name,fmt,fields):
        """Create a |BinaryParserFactory|
        
        Parameters
        ----------
        name : str
            Name for parser
        
        fmt : str
            String specifying binary format of data. See :py:mod:`struct`
        
        fields : list
            Ordered list of field names to bind to data unpacked from binary file
        """
        self.name   = name
        self.fmt    = fmt
        self.fields = fields
        self.nt = namedtuple(name,fields)
    
    def __str__(self):
        return "<%s fmt='%s' fields='%s'>" % (self.name,self.fmt,",".join(self.fields)) 
    
    def __repr__(self):
        return str(self)

    def __call__(self,fh,byte_order="<"):
        """Parse data from `fh` into a dictionary mapping field names to their values
        
        Parameters
        ----------
        fh : byte stream
            File-like pointing to binary data. Pointer in file must be
            aligned with start of record.
        
        byte_order : str
            Character indicating endian-ness of data (default: `'<'` for little-endian)
        
        Returns
        -------
        :py:class:`~collections.OrderedDict`
            Dictionary mapping field names from `self.fields` to their values
        """
        tmp_dict = self.nt._make(struct.unpack(byte_order+self.fmt,
                                               fh.read(self.calcsize(byte_order))))._asdict()
        for k in tmp_dict:
            if isinstance(tmp_dict[k],bytes):
                # Python 3.x returns bytes
                # convert byte objects to strings
                tmp_dict[k] = tmp_dict[k].decode("ascii")
            if sys.version_info < (3,) and isinstance(tmp_dict[k],unicode):
                # Python 2.x returns unicodes
                # convert unicode to strings
                tmp_dict[k] = str(tmp_dict[k].decode("ascii"))
        
        return tmp_dict

    def calcsize(self,byte_order="<"):
        """Return calculated size, in bytes, of record

        Parameters
        ----------
        byte_order : str
            Character indicating endian-ness of data (default: `'<'` for little-endian)

        Returns
        -------
        int
            Calculated size of record, in bytes
        """
        return struct.calcsize(byte_order+self.fmt)


def find_null_bytes(inp,null=b"\x00"):
    """Finds all null characters in a byte-formatted input string
    
    Parameters
    ----------
    inp : bytes (or str in Python 2.7)
        byte string
    
    Returns
    -------
    :py:class:`numpy.ndarray`
        numpy array of integers indexing where the null character *\x00* was found
    """
    indices = []
    last_found = inp.find(null)
    while last_found > -1:
        indices.append(last_found)
        last_found = inp.find(null,1+last_found)
    
    return numpy.array(indices)
