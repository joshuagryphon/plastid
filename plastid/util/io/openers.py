#!/usr/bin/env python
"""Various wrappers and utilities for opening, closing, and writing files.

Important methods
-----------------
:py:func:`argsopener`
    Opens a file for writing within a command-line script and writes to it
    all command-line arguments as a pretty-printed dictionary of metadata,
    commented out. The open file handle is then returned for subsequent
    writing

:py:func:`read_pl_table`
    Wrapper function to open a table saved by one of :data:`plastid`'s
    command-line scripts into a :class:`pandas.DataFrame`.

:py:func:`opener`
    Guesses whether a file is bzipped, gzipped, zipped, or uncompressed based upon 
    file extension, opens it appropriately, and returns a file-like object.

:py:func:`NullWriter`
    Returns an open filehandle to the system's null location.
"""
import sys
import os
import pandas as pd
from plastid.util.io.filters import AbstractWriter
from collections import Iterable

class NullWriter(AbstractWriter):
    """Writes to system-dependent null location.
    On Unix-like systems & OSX, this is typically /dev/null. On Windows, simply "nul"
    """
    
    def __init__(self):
        self.stream = open(os.devnull,"w")
    
    def filter(self,stream):
        return stream
    
    def __repr__(self):
        # unusual repr, but useful for documentation by Sphinx
        return "NullWriter()"
    
    def __str__(self):
        return self.__repr__()


def multiopen(inp,fn=None,args=None,kwargs=None):
    """Normalize filename/file-like/list of filename or file-like to a list of appropriate objects
    
    If not list-like, `inp` is converted to a list. Then, for each element `x` in
    `inp`, if `x` is file-like, it is yielded. Otherwise, `fn` is applied to `x`,
    and the result yielded. 
    
    
    Parameters
    ----------
    inp : str, file-like, or list-like of either of those
        Input describing file(s) to open
    
    fn : callable, optional
        Callable to apply to input to open it
    
    args : tuple, optional
        Tuple of positional arguments to pass to `fn`
        
    kwargs : keyword arguments
        Arguments to pass to `fn`
        
    Yields
    ------
    Object
        Result of applying `fn` to filename(s) in `inp` 
    """
    if fn is None:
        fn = lambda x, **z: x
    
    if args is None:
        args = ()
    
    if kwargs is None:
        kwargs = {}
        
    if isinstance(inp,str):
        out = [inp]
    elif isinstance(inp,Iterable):
        out = inp

    for obj in out:
#         if isinstance(obj,(safe_file,IOBase,StringIO.StringIO)): #cStringIO
#             yield obj
#         else:
#             yield fn(obj,**kwargs)
        if isinstance(obj,str):
            yield fn(obj,*args,**kwargs)
        else:
            yield obj


def opener(filename,mode="r",**kwargs):
    """Open a file, detecting whether it is compressed or not, based upon
    its file extension. Extensions are tested in the following order:
    
       +----------------+------------------+
       | File ends with | Presumed to be   |
       +================+==================+
       | gz             |    gzipped       |
       +----------------+------------------+
       | bz2            |    bzipped       |
       +----------------+------------------+
       | zip            |    zipped        |
       +----------------+------------------+
       | anything else  |    uncompressed  |
       +----------------+------------------+
        
    Parameters
    ----------
    filename : str
        Name of file to open
    
    mode : str
        Mode in which to open file. See Python standard
        libarary documentation on file opening modes for
        choices (e.g. "r", "a, "w" with or without "b")
    
    **kwargs
        Other parameters to pass to appropriate file opener 
    """
    if filename.endswith(".gz"):
        import gzip
        if "b" not in mode:
            mode += "b"
        call_func = gzip.GzipFile
    elif filename.endswith(".bz2"):
        import bz2
        if "b" not in mode:
            mode += "b"
        call_func = bz2.BZ2File
    elif filename.endswith(".zip"):
        import zipfile
        if "b" not in mode:
            mode += "b"
        call_func = zipfile.ZipFile
    else:
        call_func = open
    
    return call_func(filename,mode,**kwargs)

# needs unittest
def read_pl_table(filename,**kwargs):
    """Open a table saved by one of :data:`plastid`'s command-line scripts,
    passing default arguments to :func:`pandas.read_table`:
    
        ==========   =======
        Key          Value
        ----------   -------
        sep          `"\t"`
        comment      `"#"`
        index_col    `None`
        header       `0`
        ==========   =======
    
    Parameters
    ----------
    filename : str
        Name of file. Can be gzipped, bzipped, or zipped.
        
    kwargs : keyword arguments
        Other keyword arguments to pass to :func:`pandas.read_table`.
        Will override defaults.
        
    Returns
    -------
    :class:`pandas.DataFrame`
        Table of results
    """
    args = { "sep"      : "\t",
             "comment"  : "#",
             "index_col"  : None,
             "header"     : 0,
        }
    args.update(kwargs)
    table = pd.read_table(filename,**args)
    return table

def get_short_name(inpt,separator=os.path.sep,terminator=""):
    """Gives the basename of a filename or module name passed as a string.
    If the string doesn't match the pattern specified by the separator
    and terminator, it is returned unchanged. 
    
    Examples
    --------
    >>> get_short_name("test")
    'test'

    >>> get_short_name("test.py",terminator=".py")
    'test'

    >>> get_short_name("/home/jdoe/test.py",terminator=".py")
    'test'

    >>> get_short_name("/home/jdoe/test.py.py",terminator=".py")
    'test.py'

    >>> get_short_name("/home/jdoe/test.py.2")
    'test.py.2'
    
    >>> get_short_name("/home/jdoe/test.py.2",terminator=".py")
    'test.py.2'

    >>> get_short_name("plastid.bin.test",separator="\.",terminator="")
    'test'
    
    Parameters
    ----------
    inpt : str
        Input
    
    terminator : str
        File terminator (default: "")
    
    Returns
    -------
    str
    """
    import re
    tlen = len(terminator)
    if inpt[-tlen:] == terminator:
        inpt = inpt[:-tlen]
        
    pat = r"([^%s]+)+$" % separator
    try:
        stmp = re.search(pat,inpt).group(1)
    except AttributeError:
        return inpt
    #if terminator[0] == ".":
    #    terminator = "\"" + terminator
        
    return stmp

def argsopener(filename,namespace,mode="w",**kwargs):
    """Open a file for writing, and write to it command-line arguments
    formatted as a pretty-printed dictionary in comment metadata.
    
    Parameters
    ----------
    filename : str
        Name of file to open. If it terminates in `'.gz'` or `'.bz2'`
        the filehandle will write to a gzipped or bzipped file
                          
    namespace : :py:class:`argparse.Namespace`
        Namespace object from argparse.ArgumentParser
    
    mode : str
        Mode of writing (`'w'` or `'wb'`)
    
    **kwargs
        Other keyword arguments to pass to file opener
    
    Returns
    -------
    open filehandle
    """
    if "w" not in mode:
        mode += "w"
    fout = opener(filename,mode,**kwargs)
    fout.write(args_to_comment(namespace))
    return fout

def args_to_comment(namespace):
    """Formats a :class:`argparse.Namespace` into a comment block
    useful for printing in headers of output files
    
    Parameters
    ----------
    namespace  : :py:class:`argparse.Namespace`
        Namespace object returned by argparse.ArgumentParser
        
    Returns
    -------
    string
    """
    import datetime
    dtmp = namespace.__dict__
    ltmp = ["## date = '%s'" % datetime.datetime.today(),
            "## execstr = '%s'" % " ".join(sys.argv)
            ]
    ltmp.append("## args = {  ")
    ltmp2 = pretty_print_dict(dtmp).split("\n")[1:-2]
    for i in range(len(ltmp2)):
        ltmp2[i] = "##" + ltmp2[i]
    sout = "\n".join(ltmp) + "\n" + ",\n".join(ltmp2) + "\n##        }\n"
    return sout

def pretty_print_dict(dtmp):
    """Pretty prints an un-nested dictionary 

    Parameters
    ----------
    dtmp : dict
    
    Returns
    -------
    str
        pretty-printed dictionary
    """
    ltmp = []
    keys = dtmp.keys()
    maxlen = 2 + max([len(K) for K in keys])
    for k,v in sorted(dtmp.items(),key=lambda x: x[0]):
        if type(v) == type(""):
            v = "'%s'" % v
        new_k = "'%s'" % k
        stmp = ("          {0:<%s} : {1}," % maxlen).format(new_k,v)
        ltmp.append(stmp)
    sout = "\n".join(ltmp)
    return "{\n%s\n}\n" % sout
