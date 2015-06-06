#!/usr/bin/env python
"""Various wrappers and utilities for opening, closing, and writing files.

Important methods
-----------------
:py:func:`argsopener`
    Opens a file for writing from within a command-line script. Prepends at the
    top of the file as comment metadata a pretty-printed dictionary of the 
    command-line arguments and the date, so that each file doesn't get divorced
    from its own metadata.

:py:func:`opener`
    Guesses whether a file is bzipped, gzipped, zipped, or uncompressed based upon 
    file extension, opens it appropriately, and returns a file-like object.

:py:func:`NullWriter`
    Returns an open filehandle to the system's null location.
"""
import datetime
import gzip
import bz2
import zipfile
import re
import os
from yeti.util.io.filters import AbstractWriter

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
    
    
def guess_opener(namestub):
    """Guesses compression of potentially compressed files, based upon base 
    filename, by testing different extensions, and seeing which file is present
    and/or open. Files are sought by extension, in the following order, and the
    first one found is opened:
    
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
    namestub : str
        Name of file, without extension
    
    Returns
    -------
    filehandle pointing to open file, if found
    """
    fin = None
    suffixes = ("",".bz2",".gz",".zip")
    for suffix in suffixes:
        try:
            fin = opener("%s%s" % (namestub,suffix))
            break
        except IOError:
            continue
    return fin

def opener(filename,mode="r",**kwargs):
    """Wrapper to open compressed or uncompressed files, based upon filename, 
    as follows:
    
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
        if "b" not in mode:
            mode += "b"
        call_func = gzip.GzipFile
    elif filename.endswith(".bz2"):
        if "b" not in mode:
            mode += "b"
        call_func = bz2.BZ2File
    elif filename.endswith(".zip"):
        if "b" not in mode:
            mode += "b"
        call_func = zipfile.ZipFile
    else:
        call_func = open
    
    return call_func(filename,mode,**kwargs)
    
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

    >>> get_short_name("yeti.bin.test",separator="\.",terminator="")
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
    pat = r"([^%s]+)%s$" % (separator,terminator)
    try:
        stmp = re.search(pat,inpt).group(1)
    except AttributeError:
        return inpt
    #if terminator[0] == ".":
    #    terminator = "\"" + terminator
        
    return re.sub("%s$" % terminator,"",stmp)

def argsopener(filename,namespace,mode="w",**kwargs):
    """Opens a file for writing, formatting arguments passed in a :py:class:`argparse.Namespace`
    object as a comment quoted in the header of the file
    
    Parameters
    ----------
    filename : str
        Name of file to open. If it terminates in ".gz" or ".bz2"
        the filehandle will write to a gzipped or bzipped file
                          
    namespace : :py:class:`argparse.Namespace`
        Namespace object from argparse.ArgumentParser
    
    mode : str
        Mode of writing ('w' or 'wb')
    
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
    """Formats a Namespace from an argparse.ArgumentParser into a comment block
    useful for printing in headers of output files
    
    Parameters
    ----------
    namespace  : :py:class:`argparse.Namespace`
        Namespace object returned by argparse.ArgumentParser
        
    Returns
    -------
    string
    """
    dtmp = namespace.__dict__
    ltmp = ["## date = '%s'" % datetime.datetime.today()]
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
    str : pretty-printed dictionary, with colons aligned
    """
    ltmp = []
    keys = dtmp.keys()
    maxlen = 2 + max([len(K) for K in keys])
    for k,v in sorted(dtmp.items(),key=lambda x: x[0]):
        if type(v) == type(""):
            v = "'%s'" % v
        new_k = "'%s'" % k
        stmp = ("          {0:<%s} : {1}" % maxlen).format(new_k,v)
        ltmp.append(stmp)
    sout = "\n".join(ltmp)
    return "{\n%s\n}\n" % sout
