#!/usr/bin/env python
"""Exceptions and warnings 

Exceptions
----------
|MalformedFileError|
    Raised when a file cannot be parsed as expected, and
    execution must halt

Warnings
--------
|ArgumentWarning|
    Warning for command-line arguments that:
    
      - are nonsenical, but recoverable
      - together might cause very slow execution
        (e.g. run would be optimized by other combinations)

|FileFormatWarning|
    Warning for slightly malformed but usable files

|DataWarning|
    Warning raised when:

      - data has unexpected attributes
      - data has nonsensical, but recoverable values for attributes
      - when values are out of the domain of a given operation,
        but skipping the operation or estimating the value is permissible
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



class ArgumentWarning(Warning):
    """Warning for nonsensical but recoverable combinations of command-line arguments,
    or arguments that risk slow program execution"""
    pass


class FileFormatWarning(Warning):
    """Warning for slightly malformed but usable files"""
    pass


class DataWarning(Warning):
    """Warning for unexpected attributes of data
    Raised when:

      - data has unexpected attributes
      - data has nonsensical, but recoverable values
      - values are out of the domain of a given operation, but execution
        can continue if the value is estimated or the operation skipped
    """
