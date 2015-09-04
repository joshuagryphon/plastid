#!/usr/bin/env python
"""Utility classes, analagous to Unix-style pipes, for filtering or processing
input or output streams, such as file objects.

These filters may be composed by wrapping one around another, to perform
multiple operations in line on, for example, a text stream, before finally
coverting it to another type of data.

Various readers are included:

    - :py:class:`AbstractReader`
        the base class for all reader functions. To create a functional subclass,
        you should only need to override the :py:meth:`~AbstractReader.filter` method,
        which is called on each unit of data (usually a line of an input file).
    
    - :py:class:`CommentReader`
        read through text documents, skipping over commented lines transparently
    
    - :py:class:`SkipBlankReader`
        read through text documents, skipping over blank lines transparently
    
    - :py:class:`FunctionReader`
        apply an arbitrary function to each unit of input data
    
    - :py:class:`TeeReader`
        similar to Unix/Linux shell command ``tee``. Allows input streams to be
        sent to an arbitrary number of listeners

And various writers:

    - :py:class:`AbstractWriter`
        the base class for all reader functions. To create a functional subclass,
        you should only need to override the :py:meth:`~AbstractWriter.filter` method,
        which is called on each unit (usually a line or block of text)
    
    - :py:class:`NameDateWriter`
        Prepend timestamps to each line of string input before writing

    - :py:class:`CommentWriter`
        Filter out commented lines from text stream before writing


Examples
--------
Open a file, skipping comments and blank lines::

    >>> my_reader = CommentReader(SkipBlankReader(open("some_file.txt")))
    >>> for line in my_reader:
    >>>     pass # do something with each line, now that comments are removed

Open a file, applying a function ``foo_func`` to each line of input (*note,*
``foo_func`` can return any data type, not just a string)::

    >>> my_reader = FunctionReader(open("some_file.txt"),foo_func)
    >>> for data_unit in my_reader:
    >>>     pass # do something with each `foo`ed data unit

Write to stdout, prepending name and date::

    >>> import sys
    >>> my_writer = NameDateWriter(stream=sys.stdout)
    >>> my_writer.write(some_text)

"""
from __future__ import print_function
import sys
import datetime
from abc import abstractmethod

#===============================================================================
# INDEX: readers
#===============================================================================

class AbstractReader(object):
    """Abstract base class for stream-reading filters. These may be wrapped around
    open-file like objects, for example, to remove comments or blank lines from text,
    or to convert units of input from one data type to another.
    
    Create a filter by subclassing this, and defining `self.filter()`
    
    See also
    --------
    CommentReader
        A reader that removes comments from text data
    
    """
    
    def __init__(self,stream):
        """Create an |AbstractReader|
        
        Parameters
        ----------
        stream : file-like
            Input data
        """
        self.stream=stream
    
    def __del__(self):
        try:
            self.close()
        except AttributeError:
            pass
        
    def __next__(self):
        return self.filter(next(self.stream))

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def read(self):
        """Similar to :py:func:`file.read`. Process all units of data, assuming it is string-like
        
        Returns
        -------
        str
        """
        return "".join(self.readlines())

    def readline(self):
        """Similar to :py:func:`file.readline`. You probably want to use next(self) instead.
        Process a single line of data, assuming it is string-like
        
        Returns
        -------
        object
            a unit of processed data
        """
        return self.filter(self.stream.readline())

    def readlines(self):
        """Similar to :py:func:`file.readlines`.
        
        Returns
        -------
        list
            processed data
        """
        lreturn = []
        for line in self:
            lreturn.append(line)

        return lreturn
    
    def close(self):
        """Close stream"""
        self.stream.close()

    @abstractmethod
    def filter(self,data):
        """Method that filters or processes each unit of data.
        Override this in your subclass
        
        Parameters
        ----------
        data : unit of data
            Whatever data to filter/format. Often string, but not necessary
        
        Returns
        -------
        object
            formatted data. Often string, but not necessarily
        """
        pass


class FunctionReader(AbstractReader):
    """Apply a function to each line in an input stream"""
    
    def __init__(self,stream,func):
        """Create a FunctionReader
        
        Parameters
        ----------
        stream : file-like
            Input stream
            
        func : function
            Function to apply to each unit of input in `stream`
        """
        self.filter = func
        AbstractReader.__init__(self,stream)


class SkipBlankReader(AbstractReader):
    """Ignores blank/whitespace-only lines in a text stream"""

    def filter(self,line):
        """Return next non-blank line of text
        
        Parameters
        ----------
        line : str
            Line of text

        Returns
        -------
        str
        """
        if len(line.strip()) == 0:
            return self.__next__()
        else:
            return line
        
        
class CommentReader(AbstractReader):
    """Ignore lines beginning with `'#'`, optionally preceded by whitespace
    from a text stream. Comments beginning mid line are left in-place
    """

    def __init__(self,stream):
        """Create a CommentReader
        
        Parameters
        ----------
        stream : file-like
            Input data
        """
        self.comments = []
        AbstractReader.__init__(self,stream)

    def get_comments(self):
        """Return all of the comments that have been found so far.
        Does not remove comments from the internal list.

        Returns
        -------
        list
            Comments found in text
        """
        return self.comments
    
    def filter(self,line):
        """Return next non-commented line of text
        
        Parameters
        ----------
        line : str
            Line of text

        Returns
        -------
        str
        """
        ltmp = line.lstrip()
        if len(ltmp) > 1 and ltmp[0] == "#":
            self.comments.append(line.strip())
            return self.__next__()
        else:
            return line


class BackwardReader(AbstractReader):
    """Reverses each line of a text stream character-wise."""
    def filter(self,line):
        """Return next non-commented line of text
        
        Parameters
        ----------
        line : str
            Line of text

        Returns
        -------
        str
            Reversed line of text
        """        
        return line [::-1]


class TeeReader(AbstractReader):
    """Tee an input stream to listeners that register themselves with the TeeReader
    via the add_listener() method. Similar to shell command ``tee``

    Each listener must implement an ``alert()`` method, in order to receive the data.
    If `alert()` is not implemented, errors are suppressed.
    
    See also
    --------
    TeeListener : an example of a listener class
    """

    def __init__(self,stream):
        """Create an TeeReader
        
        Parameters
        ----------
        stream : file-like
            Input data
        """
        self.listeners = []
        AbstractReader.__init__(self,stream)

    def add_listener(self,listener):
        """Register a single listener with this reader
        
        Parameters
        ----------
        listener : |TeeListener|-like
        """
        self.listeners.append(listener)
    
    def add_listeners(self,*many_listeners):
        """Register one or more listeners with this reader
        
        Parameters
        ----------
        many_listeners : one or more |TeeListener|-like
        """
        for listener in many_listeners:
            self.add_listener(listener)
    
    def filter(self,line):
        """Sends each line to each listener. Complains if listener cannot listen!
        
        Parameters
        ----------
        line : a unit of input data
        
        Returns
        -------
        input data
        """
        for listener in self.listeners:
            try:
                listener.alert(line)
            except AttributeError:
                import warnings
                warnings.warn("Could not alert listener %s: " % str(listener))
        return line
    
    
class TeeListener(object):
    """Listener class for TeeFilter. Listeners if registered with a |TeeReader|
    will receive and process each unit of input via its ``alert()`` method"""
    
    @abstractmethod
    def alert(self,data):
        """Process input from a |TeeReader|.
        Override this method to implement whatever function you need.
        
        Parameters
        ----------
        data : object
            A unit of data
        """
        pass

class TestTeeListener(TeeListener):
    """Example of a TeeListener"""
    
    def alert(self,line):
        print(self.name + " heard something: "+line)


#===============================================================================
# INDEX: writers
#===============================================================================

class AbstractWriter(object):
    """Abstract base class for stream-writing filters.
    Create a filter by subclassing this, and defining self.filter().
    
    See also
    --------
    NameDateWriter : a writer that prepends timestamps to text
    """
    def __init__(self,stream):
        """Create an AbstractWriter
        
        Parameters
        ----------
        stream : file-like, open for writing
            Output stream to which filtered/formatted data will be written
        """
        self.stream=stream
    
    def __del__(self):
        """Safely close streams and then destroy |AbstractWriter|
        """
        try:
            self.close()
        except:
            pass
            
    def write(self,data):
        """Write data to `self.stream`
        
        Parameters
        ----------
        data : unit of data
            Whatever data to filter/format. Often string, but not necessary
        """
        self.stream.write(self.filter(data))
    
    def flush(self):
        """Flush `self.stream`"""
        self.stream.flush()

    def close(self):
        """flush and close `self.stream`"""
        self.flush()
        self.stream.close()
    
    @abstractmethod
    def filter(self,data):
        """Method that filters or processes each unit of data.
        Override this in your subclass
        
        Parameters
        ----------
        data : unit of data
            Whatever data to filter/format. Often string, but not necessarily
        
        Returns
        -------
        object
            formatted data. Often string, but not necessary
        """ 
        pass


class NameDateWriter(AbstractWriter):
    """Prepend date and program name to each line of output"""
    def __init__(self,name,line_delimiter="\n",stream=sys.stderr):
        """Create a NameDateWriter
        
        Parameters
        ----------
        name : str
            Name to prepend
        
        stream : file-like
            Stream to write to (Default: :obj:`sys.stderr`)
        
        line_delimiter : str, optional
            Delimiter, postpended to lines. (Default `'\n'`)
        """
        self.name = name
        self.delimiter = line_delimiter
        AbstractWriter.__init__(self,stream)
    
    def filter(self,line):
        """Prepend date and time to each line of input
        
        Parameters
        ----------
        line : str
            Input
        
        Returns
        -------
        str : Input with date and time prepended
        """
        now = datetime.datetime.now()
        return "%s [%s]: %s%s" % (self.name, now, line.strip(self.delimiter), self.delimiter)
    
    # included for backward compatibility
    def __call__(self,line):
        self.write(line)


class CommentWriter(AbstractWriter):
    """Filter out lines beginning with `'#'` from data written to a stream"""

    def filter(self,line):
        if line.startswith("#"):
            return ""
        else:
            return line

# alias included for backward compatibility
Printer = NameDateWriter
