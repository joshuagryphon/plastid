#!/usr/bin/env python
"""Exceptions and warnings used by command-line scripts or other methods

Exceptions
----------
|MalformedFileError|
    Raised when a file cannot be parsed as expected

"""

class MalformedFileError(Exception):
    """Exception class for when files cannot be parsed as they should be
    """
    
    def __init__(self,filename,message,line_num=None):
        """Create a |MalformedFileError|
        
        Parameters
        ----------
        filename : str
            Name of file causing problem
        
        message : str
            Message explaining how the file is malformed.
        
        line_num : int or None, optional
            Number of line causing problems
        """
        self.filename = filename
        self.msg      = message
        self.line_num = line_num
    
    def __str__(self):
        if self.line is None:
            return "Error opening file '%s': %s" % (self.filename, self.msg)
        else:
            return "Error opening file '%s' at line %s: %s" % (self.filename, self.line_num, self.msg)


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


