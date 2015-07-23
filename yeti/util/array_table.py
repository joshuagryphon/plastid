#!/usr/bin/env python
"""Describes an |ArrayTable| class, which is essentially a big spreadsheet
build on top of :class:`numpy.ndarray` objects.

Notes
-----
This class is deprecated in favor of :class:`pandas.DataFrame` and will gradually
be replaced by it.
"""
import copy
import numpy
import pandas as pd
from yeti.util.services.misc import guess_formatter
from yeti.util.io.openers import NullWriter

_NUMERIC_DTYPES = "iufc"


class ArrayTable(object):
    """A representation of a spreadsheet, in which each column has a single
    data type. Multiple data types are allowed within single sheets. Internally,
    the spreadsheet is represented as a dictionary of :py:class:`numpy.ndarray` s. For
    convenience, access to this dictionary is exposed. Note: the length of
    the |ArrayTable| is defined as the number of rows in the table, NOT the number
    of columns (as is the case with a dictionary).
    
    Single columns of the spreadsheet my be accessed using the corresponding
    header key. Rows may be accessed using get_rows() in combination with a 
    selection mask.
    
    While this class can be instantiated directly from a dictionary of
    :py:class:`numpy.ndarray` objects, it is more convenient to use it to parse
    tab-delimited files using :meth:`ArrayTable.from_file`,
    and to write such files using :meth:`ArrayTable.to_file`
    """
    def __init__(self,my_dict):
        """Create an |ArrayTable| from a dictionary of sequences
        
        Parameters
        ----------
        my_dict : dict
            Dictionary of sequences, each corresponding to a column
            in the |ArrayTable|
        
        Raises
        ------
        |ValueError| if sequences have different lengths
        """
        lengths = set([len(X) for X in my_dict.values()])
        if len(lengths) != 1:
            length_list = ["%s\t%s" % (K,len(V)) for K, V in my_dict.items() ]
            raise ValueError("Lengths of all columns must be identical: %s" % "\n".join([X for X in length_list]))
        # convert to array in case incoming dict was lists or, worse, MaskedArrays
        # (MaskedArrays, when written to file, convert nans to 0!!)
        self._dict = { K : numpy.array(V) for K,V in my_dict.items() }
    
    def keys(self):
        """Return names of columns in |ArrayTable|"""
        return self._dict.keys()
    
    def values(self):
        """Return sequence of columns in |ArrayTable|"""
        return self._dict.values()
    
    def items(self):
        """Return sequence of (key,column) pairs for column names and columns
        in |ArrayTable|
        """
        return self._dict.items()
    
    def shape(self):
        """Return a tuple of number of rows and colums in the dictionary"""
        return (len(self),len(self.keys()))
    
    def __eq__(self,other):
        """Test equality with another |ArrayTable|, defined as equality of keys,
        equality of lengths, identity of all non-numerical values, and equality
        of all numerical values within a tolerance
        
        Parameters
        ----------
        other : |ArrayTable|
        
        Returns
        -------
        bool
        
        See Also
        --------
        equal_enough
            column-wise equality metric, allowing for slight differences in floats
        """
        if len(self) != len(other):
            return False
        elif set(self.keys()) != set(other.keys()):
            return False
        else:
            match = True
            for k in self.keys():
                match = match & equal_enough(self[k],other[k])
                if match is False:
                    return match
        return match
        
    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "<ArrayTable %sx%s>" % self.shape()
    
    def __len__(self):
        """Return the number of rows in the |ArrayTable|
        
        Returns
        -------
        int
            number of rows in the |ArrayTable|
        """
        return len(list(self.values())[0])

    def __iter__(self):
        """Enable iteration- defined over keys as for a `dict`"""
        return iter(sorted(list(self.keys())))

    def __getitem__(self,key):
        """Retrieve a single column from the |ArrayTable|
        
        Parameters
        ----------
        key : str
            Key corresponding to desired column
        
        Returns
        -------
        :py:class:`numpy.ndarray`
            numpy array representing a column of data
        """
        return self._dict[key]
    
    def __setitem__(self,key,val):
        """Set a single column in the |ArrayTable|
        
        Parameters
        ----------
        key : str
            Key corresponding to new column
        
        val : array-like
            New column, as a sequence
        
        Raises
        ------
        ValueError
            if `len(val) != len(self)`
        """
        if len(val) != len(self):
            raise ValueError("Length of incoming array (%s) is unequal length of destination table (%s)" % (len(val),len(self)))
        self._dict[key] = numpy.array(val)
    
    def get_rows(self,keyorder=None,mask=None):
        """Return a section of the |ArrayTable| as a set of rows
        
        Parameters
        ----------
        keyorder : list
            Sequence of names corresponding to columns from which to draw values
                     
        mask : array-like
            logical mask or sequence of indices specifying rows from which to draw values
        
        Returns
        -------
        list<tuple> of values
        """
        if keyorder is None:
            keyorder = self.keys()
        if mask is None:
            mask = numpy.tile(True,len(self))
        return zip(*tuple([self[K][mask] for K in keyorder]))
    
    def merge(self,source,host_match_key,source_match_key,keys_to_copy=None):
        """Merge columns from `source` into an existing |ArrayTable| in-place.
        `source` is not required to have the same length as `SELF`.
        Missing data are treated as `numpy.nan` for numerical types, or `None` 
        for non-numerical types.
        
        Parameters
        ----------
        source : |ArrayTable|
            |ArrayTable| from which data will be copied
            
        host_match_key : str
            Key in host |ArrayTable| on which to perform join
            
        source_match_key : str
            Key in source |ArrayTable| on which to perform join
            
        keys_to_copy : list<str> or `None`
            Sequence of keys indicating which columns to copy
            from source to host. If `None`, all keys in source,
            save the match key, are copied.
        """
        if keys_to_copy is None:
            keys_to_copy = set(source.keys()) - { source_match_key }
        
        my_len = len(list(self.values())[0])
        for k in keys_to_copy:
            if source[k].dtype.kind in _NUMERIC_DTYPES:
                self[k] = numpy.tile(numpy.nan,my_len) 
            else:
                self[k] = numpy.tile(None,my_len) 
            
        for i in range(len(list(source.values())[0])):
            match_mask = (self[host_match_key] == source[source_match_key][i])
            for k in keys_to_copy:
                self[k][match_mask] = source[k][i]
    
    def as_array(self,keyorder=None,mask=None):
        """Return a portion of the |ArrayTable| as an MxN ``numpy`` array.
        Note, all keys specified in keyorder must be numerical types. Otherwise,
        ``numpy`` will throw an error.
        
        Parameters
        ----------
        keyorder : list<str>, optional
            sequence of keys reflecting desired column order in output file.
            If `None`, keys are sorted alphabetically
                         
        mask : array-like, optionals
            logical mask or sequence of indices specifying rows from which to draw values

        Returns
        -------
        :py:class:`numpy.ndarray`
        """
        return numpy.array(list(self.get_rows(keyorder=keyorder,mask=mask)))
    
    def to_file(self,fh,keyorder=None,mask=None,formatters={}):
        """Writes a dictionary of numpy arrays to a tab-delimited text file.
        For floats, forces 64-bit precision. Integars are explicitly recorded
        as integers. All other datatypes are written as str
        
        Parameters
        ----------
        fh : file-like
            filehandle to write to
        
        keyorder : list<str>, optional
            sequence of keys reflecting desired column order in
            output file. if `None`, keys are sorted alphabetically
                         
        mask : array-like, optional
            logical mask or sequence of indices specifying rows from which to draw values
        
        formatters : dict, optional
            Override default formatting of an object. Supply the object class
            as a key, and a string formatter as a value. For example `str`
            or `'{:d}'.format`
        """

        format_dict = { numpy.float128  : '{:.128e}'.format,
                        numpy.float64   : '{:.64e}'.format,
                        numpy.float32   : '{:.32e}'.format,
                        numpy.float16   : '{:.16e}'.format,
                        numpy.int64     : '{:d}'.format,
                        numpy.int32     : '{:d}'.format,
                        numpy.uint      : '{:d}'.format,
                        numpy.uint0     : '{:d}'.format,
                        numpy.uint8     : '{:d}'.format,
                        numpy.uint16    : '{:d}'.format,
                        numpy.uint32    : '{:d}'.format,
                        "default"       : str,
                        object          : str,
                        
                       }
        
        for k,v in formatters.items():
            format_dict[k] = v
        
        if keyorder is None:
            keyorder = sorted(self.keys())
            
        fh.write("#" + "\t".join(keyorder) + "\n")
        dtypes = [self[K].dtype.type for K in keyorder]
        formatters = [format_dict.get(X,format_dict["default"]) for X in dtypes]
        for row in self.get_rows(keyorder=keyorder,mask=mask):
            ltmp = [F(X) for F,X in zip(formatters,row)]
            fh.write("\t".join(ltmp) + "\n")
    
    def as_data_frame(self):
        """Transitional convenience method to convert to a :py:class:`pandas.DataFrame`
        
        Returns
        -------
        :py:class:`pandas.DataFrame`
        """
        return pd.DataFrame(self._dict)
    
    @staticmethod
    def from_file_to_data_frame(fh):
        """Transitional convenience method to read a file from :py:meth:`ArrayTable.to_file`
        into a :class:`pandas.DataFrame`
        
        Parameters
        ----------
        fh : file-like
            filehandle pointing to data
        
        Returns
        -------
        :py:class:`pandas.DataFrame`
        """
        # inefficient but functional.
        return ArrayTable.from_file(fh).as_data_frame()
        
    @staticmethod
    def concat(*arraytables):
        """Concatenate multiple |ArrayTables| into a single |ArrayTable| in order.
        
        Parameters
        ----------
        arraytables : one or more |ArrayTables|
            sequence of |ArrayTables| to merge
            
        Returns
        -------
        |ArrayTable|
            merged |ArrayTable|
            
        Raises
        ------
        ValueError : if input |ArrayTables| have mismatched columns
        """
        keys = set()
        for at in arraytables:
            keys |= set(at.keys())
        for at in arraytables:
            if set(at.keys()) != keys:
                raise ValueError("ArrayTable columns don't match!")

        dtmp = {}
        for k in keys:
            dtmp[k] = []
            for at in arraytables:
                dtmp[k].extend(list(at[k]))
        
        dtmp = { K : numpy.array(V) for K,V in dtmp.items() }
        dtmp = ArrayTable(dtmp)
        return dtmp
    
    @staticmethod
    def from_file(fh,formatters=guess_formatter):
        """Read a tab-delimited file into an |ArrayTable|.
        The first line that begins with a single pound sign is taken to be a
        header row; these headers will be the keys in the returned dictionary.
        Any line beginning with two pound signs (`'##'`) is assumed to contain
        metadata, which is discarded. 
        
        Parameters
        ----------
        fh : file-like
            filehandle pointing to data
        
        formatters
            a sequence of formatters (e.g. `str`, `int`, `float`) used to parse
            each column, or a single formatter to apply to all columns
            (Default: guess which formatter to use)
                            
        Returns
        -------
        dict[key] = numpy.array
        """
        headers = []
        dtmp = {}
        headers_done = False
        for line in fh:
            if line.startswith("##"):
                pass
            elif line.startswith("#"):
                if headers_done is True:
                    continue
                else:
                    headers = line.strip("\n")[1:].split("\t")
                    dtmp = { K : copy.deepcopy([]) for K in headers }
                    headers_done = True
                    if type(formatters) not in (type([]),type(())):
                        formatters = [formatters] * len(headers)
                    else:
                        formatters = formatters
            else:
                items = line.strip("\n").split("\t")
                for k,v,formatter in zip(headers,items,formatters):
                    dtmp[k].append(formatter(v))
        dtmp = { K : numpy.array(V) for K,V in dtmp.items() }
        return ArrayTable(dtmp)


def equal_enough(col1,col2,tol=1e-10,printer=NullWriter()):
        """If `col1` and `col2` are both numeric, test that
        all their values are within `tol` of each other. :obj:`numpy.nan` values,
        if present, must be in the same place in each column. Ditto :obj:`numpy.inf`
        values.
        
        If `col1` and `col2` are not numeric, return true if they have the same
        `dtype` and the same values in all cells.
        
        Parameters
        ----------
        col1 : numpy.ndarray
            First column of data
        
        col2 : numpy.ndarray
            Second column of data
            
        tol : float
            Error tolerance for numeric data between col1 and col2, value-wise
            
        printer : anything implementing a ``write()`` method (e.g. a |NameDateWriter|)
            if not None, rich comparison information will be sent to this writer
        
        Returns
        -------
        bool
            `True` if `col1 == col2` for non-numeric data;
            `True` if `abs(col1 - col2) <= tol` for numeric data;
            `False` otherwise
        """
        dtype_test  = col1.dtype.kind in _NUMERIC_DTYPES
        dtype_test &= col2.dtype.kind in _NUMERIC_DTYPES
        
        if dtype_test:
            nan1 = numpy.isnan(col1)
            nan2 = numpy.isnan(col2)
            
            inf1 = numpy.isinf(col1)
            inf2 = numpy.isinf(col2)
            
            # check nans in same place
            nan_test = (nan1 == nan2).all()

            # make sure infs are in same place, and are same sign
            inf_test = (inf1 == inf2).all() & (col1[inf1] == col2[inf1]).all()
        
            diff_test = (abs(col1[~nan1 & ~inf1] - col2[~nan1 & ~inf1]) <= tol).all()
            
            if not nan_test:
                printer.write("Failed nan test")
            if not inf_test:
                printer.write("Failed inf test")
            if not diff_test:
                printer.write("Failed tolerance test")
        
            return nan_test & diff_test & inf_test
        elif col1.dtype.kind not in _NUMERIC_DTYPES and col2.dtype.kind not in _NUMERIC_DTYPES:
            diff_test = (col1 == col2)
            if not diff_test.all():
                ltmp = ["%s,%s,%s" % (N,val1,val2) for N,(val1,val2) in enumerate(zip(col1[~diff_test],col2[~diff_test]))] 
                printer.write("Differences:\n%s" % "\n".join(ltmp))
            return diff_test.all()
        else:
            return False
