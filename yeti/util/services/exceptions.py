#!/usr/bin/env python
"""Exceptions used by command-line scripts or other methods

Exceptions
----------
|MalformedFileError|
    Raised when a file cannot be parsed as expected
"""

class MalformedFileError(Exception):
    """Exception class for when files cannot be parsed as they should be
    """
    
    def __init__(self,filename,message,line_num=None):
        """Create a |MalformedFileException|
        
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
